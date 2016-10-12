#include "vtk_output.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include <stdlib.h>
#include "mkl_cblas.h"

#include "PGFEM_io.h"
#include "enumerations.h"
#include "gen_path.h"
#include "utils.h"
#include "allocation.h"
#include "interface_macro.h"
#include "PGFem3D_to_VTK.hpp"

static const char *out_dir = "VTK";
static const char *step_dir = "STEP_";
static const int ndim = 3;

void VTK_print_master(char *path,
          char *base_name,
          int time,
          int nproc,
          const PGFem3D_opt *opts)
{
  /* Print the master file wich points to all the data files and
     declares the datatypes */
  char cur_dir = '.';
  char *ptr_path = path;
  char dir_name[500];

  const int analysis = opts->analysis_type;

  if(ptr_path == NULL) ptr_path = &cur_dir;
  sprintf(dir_name,"%s/%s",ptr_path,out_dir);

  if(make_path(dir_name,DIR_MODE) != 0){
    PGFEM_printf("Directory (%s) not created!\n",dir_name);
    abort();
  }

  /* Build filename and open file */
  char *filename;
  filename = (char*) PGFEM_calloc(500, sizeof(char));
  sprintf(filename,"%s/%s_%d.pvtu",dir_name,base_name,time);
  FILE *out;
  if((out = fopen(filename,"w")) == NULL){
    PGFEM_printerr("ERROR: Cannot open file to write pvtu file! (%s)\n",
      filename);
    abort();
  }

  sprintf(dir_name,"%s/%s/%s%.5d",ptr_path,out_dir,step_dir,time);

  /* write header information */
  PGFEM_fprintf(out,"<?xml version=\"1.0\"?>\n");
  PGFEM_fprintf(out,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\""
    " byte_order=\"LittleEndian\">\n");
  PGFEM_fprintf(out,"<PUnstructuredGrid>\n");

  /* write datatypes */
  /* Point Data */
  PGFEM_fprintf(out,"<PPointData>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"Displacement\""
    " NumberOfComponents=\"3\"/>\n");

  if(opts->multi_scale){
    PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"MacroDisplacement\""
      " NumberOfComponents=\"3\"/>\n");
  }

  if(analysis == STABILIZED
     || analysis == MINI
     || analysis == MINI_3F){
    PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"Pressure\"/>\n");
  }
  PGFEM_fprintf(out,"</PPointData>\n");

  /* Cell Data */
  PGFEM_fprintf(out,"<PCellData>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"CauchyStress\""
    " NumberOfComponents=\"6\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"EulerStrain\""
    " NumberOfComponents=\"6\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"EffectiveStrain\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"EffectiveStress\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Int64\" Name=\"CellProperty\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"Damage\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"Chi\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"F\""
    " NumberOfComponents=\"9\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"P\""
    " NumberOfComponents=\"9\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"W\"/>\n");
  PGFEM_fprintf(out,"</PCellData>\n");

  /* Points */
  PGFEM_fprintf(out,"<PPoints>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
  PGFEM_fprintf(out,"</PPoints>\n");

  /* write piece list */
  sprintf(dir_name,"%s%.5d",step_dir,time);
  for (int i=0; i< nproc; i++){
    sprintf(filename,"%s/%s_%d_%d.vtu",dir_name,base_name,i,time);
    PGFEM_fprintf(out,"<Piece Source=\"%s\"/>\n",filename);
  }
  /* write footer */
  PGFEM_fprintf(out,"</PUnstructuredGrid>\n");
  PGFEM_fprintf(out,"</VTKFile>\n");

  /* close file */
  fclose(out);
  free(filename);
}

void VTK_print_cohesive_master(char *path,
             char *base_name,
             int time,
             int nproc,
             const PGFem3D_opt *opts)
{
  /* Print the master file wich points to all the data files and
     declares the datatypes */

  char cur_dir = '.';
  char *ptr_path = path;
  char dir_name[500];

  if(ptr_path == NULL) ptr_path = &cur_dir;
  sprintf(dir_name,"%s/%s",ptr_path,out_dir);

  if(make_path(dir_name,DIR_MODE) != 0){
    PGFEM_printf("Directory (%s) not created!\n",dir_name);
    abort();
  }    

  /* Build filename and open file */
  char *filename;
  filename = (char*) PGFEM_calloc(500, sizeof(char));
  sprintf(filename,"%s/%s_cohesive_%d.pvtu",dir_name,base_name,time);
  FILE *out;
  if((out = fopen(filename,"w")) == NULL){
    PGFEM_printerr("ERROR: Cannot open file to write pvtu file! (%s)\n",
      filename);
    abort();
  }

  sprintf(dir_name,"%s/%s/%s%.5d",ptr_path,out_dir,step_dir,time);

  /* write header information */
  PGFEM_fprintf(out,"<?xml version=\"1.0\"?>\n");
  PGFEM_fprintf(out,"<VTKFile type=\"PUnstructuredGrid\" "
    "version=\"0.1\" byte_order=\"LittleEndian\">\n");
  PGFEM_fprintf(out,"<PUnstructuredGrid>\n");

  /* write datatypes */
  /* Point Data */
  PGFEM_fprintf(out,"<PPointData>\n");
  if(opts->multi_scale){
    PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"MacroDisplacement\""
      " NumberOfComponents=\"3\"/>\n");
  }
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"Displacement\""
    " NumberOfComponents=\"3\"/>\n");
  PGFEM_fprintf(out,"</PPointData>\n");

  /* Cell Data */
  PGFEM_fprintf(out,"<PCellData>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Int64\" Name=\"Property\""
    " NumberOfComponents=\"1\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"OpeningDisp\""
    " NumberOfComponents=\"3\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"Tractions\""
    " NumberOfComponents=\"3\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"EffOpening\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"NormalOpening\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"ShearOpening\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"EffTraction\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"NormalTraction\"/>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float64\" Name=\"ShearTraction\"/>\n");
  PGFEM_fprintf(out,"</PCellData>\n");

  /* Points */
  PGFEM_fprintf(out,"<PPoints>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
  PGFEM_fprintf(out,"</PPoints>\n");

  /* write piece list */
  sprintf(dir_name,"%s%.5d",step_dir,time);
  for (int i=0; i< nproc; i++){
    sprintf(filename,"%s/%s_cohesive_%d_%d.vtu",dir_name,base_name,i,time);
    PGFEM_fprintf(out,"<Piece Source=\"%s\"/>\n",filename);
  }
  /* write footer */
  PGFEM_fprintf(out,"</PUnstructuredGrid>\n");
  PGFEM_fprintf(out,"</VTKFile>\n");

  /* close file */
  fclose(out);
  free(filename);
}

/****  CURRENTLY ONLY SUPPORT STB ELEMENTS  ****/
void VTK_print_vtu(char *path,
       char *base_name,
       int time,
       int myrank,
       long ne,
       long nn,
       NODE *node,
       ELEMENT *elem,
       SUPP sup,
       double *r,
       SIG *sig,
       EPS *eps,
       const PGFem3D_opt *opts,
       const int mp_id)
{
  char cur_dir = '.';
  char *ptr_path = path;
  char dir_name[500];
  if(ptr_path == NULL) ptr_path = &cur_dir;
  sprintf(dir_name,"%s/%s/%s%.5d",ptr_path,out_dir,step_dir,time);
  if(make_path(dir_name,DIR_MODE) != 0){
    PGFEM_printf("Directory (%s) not created!\n",dir_name);
    abort();
  }

  /* Build filename and open file */
  char *filename;
  filename = (char*) PGFEM_calloc(500, sizeof(char));
  sprintf(filename,"%s/%s_%d_%d.vtu",dir_name,base_name,myrank,time);

#ifndef NO_VTK_LIB
  /* USE VTK LIBRARY FOR OUTPUT */
  void *grid = PGFem3D_to_vtkUnstructuredGrid(nn,ne,node,elem,sup,sig,
                eps,r,opts->analysis_type);
  /* ASCII */
  /* PGFem3D_write_vtkUnstructuredGrid(filename,grid,1); */
  /* BINARY */
  PGFem3D_write_vtkUnstructuredGrid(filename,grid,opts->ascii);

  PGFem3D_destroy_vtkUnstructuredGrid(grid);
  free(filename);
#else

  FILE *out;
  if((out = fopen(filename,"w")) == NULL){
    PGFEM_printerr("ERROR: Cannot open file to write pvtu file! (%s)\n",
      filename);
    abort();
  }
  free(filename);

  /* Print header information */
  PGFEM_fprintf(out,"<?xml version=\"1.0\"?>\n");
  PGFEM_fprintf(out,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
    " byte_order=\"LittleEndian\">\n");
  PGFEM_fprintf(out,"<UnstructuredGrid>\n");
  PGFEM_fprintf(out,"<Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\" >\n",nn,ne);

  /* Point Data */
  PGFEM_fprintf(out,"<PointData>\n");

  /* Displacement */
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"Displacement\""
    " NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (int i=0; i<nn; i++){ /* For each node */
    for (int j=0; j<3; j++){ /* For each direction */
      if (node[i].id_map[mp_id].id[j] == 0) PGFEM_fprintf(out,"%12.12e ",0.0);
      if (node[i].id_map[mp_id].id[j] >  0) PGFEM_fprintf(out,"%12.12e ",r[node[i].id_map[mp_id].id[j]-1]);
      if (node[i].id_map[mp_id].id[j] <  0) PGFEM_fprintf(out,"%12.12e ",sup->defl[abs(node[i].id_map[mp_id].id[j])-1]);
    }
    PGFEM_fprintf(out,"\n");
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  /* MacroDisplacement */
  if(sup->multi_scale){
    PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"MacroDisplacement\""
      " NumberOfComponents=\"3\" format=\"ascii\">\n");
    double *jump = aloc1(3);
    compute_macro_grad_u(sup->F0,sup,opts->analysis_type);
    for(int i=0; i<nn; i++){
      compute_interface_macro_disp_at_node(jump,&node[i],sup->F0,opts->analysis_type);
      PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",jump[0],jump[1],jump[2]);
    }
    free(jump);
    PGFEM_fprintf(out,"</DataArray>\n");
  }

  /* Pressure */
  switch(opts->analysis_type){
  case TF:
    if(elem[0].toe==10)
      break;
  case STABILIZED:
  case MINI:
  case MINI_3F:
    PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n");
    for (int i=0; i<nn; i++){
      PGFEM_fprintf(out,"%12.12e\n",r[node[i].id_map[mp_id].id[3]-1]);
    }
    PGFEM_fprintf(out,"</DataArray>\n");
    break;
  default:
    break;
  }

  /* Cell Data */
  PGFEM_fprintf(out,"</PointData>\n");
  PGFEM_fprintf(out,"<CellData>\n");

  /* Stress */
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"CauchyStress\""
    " NumberOfComponents=\"6\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e\n",
      sig[i].el.o[0],sig[i].el.o[1],sig[i].el.o[2],
      sig[i].el.o[5],sig[i].el.o[3],sig[i].el.o[4]);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Strain */
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"EulerStrain\""
    " NumberOfComponents=\"6\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e\n",
      eps[i].el.o[0],eps[i].el.o[1],eps[i].el.o[2],
      eps[i].el.o[5]/2,eps[i].el.o[3]/2,eps[i].el.o[4]/2);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Effective Strain */
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"EffectiveStrain\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%12.12e\n",eps[i].el.eq);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Effective Stress */
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"EffectiveStress\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%12.12e\n",sig[i].el.eq);
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  /* Material */
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"CellProperty\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%ld\n",elem[i].pr);
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  
  /* Damage parameter(s) */
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" "
    "Name=\"Damage\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%12.12e\n",eps[i].dam[0].wn); /* Linear element has 1 ip */
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"Chi\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    PGFEM_fprintf(out,"%12.12e\n",eps[i].dam[0].Xn); /* Linear element has 1 ip */
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"F\""
    " NumberOfComponents=\"9\" format=\"ascii\">\n");
  for(int i = 0; i < ne; i++){
    const double *F = eps[i].il[0].F;
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e"
      " %12.12e %12.12e %12.12e %12.12e %12.12e\n",
      *(F+0),*(F+1),*(F+2),*(F+3),*(F+4),*(F+5),
      *(F+6),*(F+7),*(F+8));
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"P\""
    " NumberOfComponents=\"9\" format=\"ascii\">\n");
  {
    double *S = PGFEM_calloc(9,sizeof(double));
    double *P = PGFEM_calloc(9,sizeof(double));
    for(int i = 0; i < ne; i++){
      const double *F = eps[i].il[0].F;
      S[0] = sig[i].il[0].o[0];
      S[1] = sig[i].il[0].o[5];
      S[2] = sig[i].il[0].o[4];
      
      S[3] = sig[i].il[0].o[5];
      S[4] = sig[i].il[0].o[1];
      S[5] = sig[i].il[0].o[3];

      S[6] = sig[i].il[0].o[3];
      S[7] = sig[i].il[0].o[4];
      S[8] = sig[i].il[0].o[2];

      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ndim,ndim,ndim,
      1.0,F,ndim,S,ndim,0.0,P,ndim);

      PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e"
        " %12.12e %12.12e %12.12e %12.12e %12.12e\n",
        *(P+0),*(P+1),*(P+2),*(P+3),*(P+4),*(P+5),
        *(P+6),*(P+7),*(P+8));
    }
    free(S);
    free(P);
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  
  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"W\" format=\"ascii\">\n");
  for(int i = 0; i < ne; i++){
    PGFEM_fprintf(out,"%12.12e\n",eps[i].il[0].Y);
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  
  if(opts->analysis_type==TF && elem[0].toe==10)
  {  
    /* Pressure */
    PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"TF_Pressure\" format=\"ascii\">\n");
    for (int i=0; i<ne; i++){
      PGFEM_fprintf(out,"%12.12e\n",eps[i].d_T[0]);
    }
    PGFEM_fprintf(out,"</DataArray>\n");  
  }  
  
  if(opts->analysis_type==TF)
  {  
    /* Volume */
    PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"TF_Volume\" format=\"ascii\">\n");
    for (int i=0; i<ne; i++){
      PGFEM_fprintf(out,"%12.12e\n",eps[i].T[0]);
    }
    PGFEM_fprintf(out,"</DataArray>\n");  
  }
  

  /* End Cell data */
  PGFEM_fprintf(out,"</CellData>\n");

  /* Nodes */
  PGFEM_fprintf(out,"<Points>\n");
  PGFEM_fprintf(out,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0; i<nn; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",
      node[i].x1_fd,node[i].x2_fd,node[i].x3_fd);
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  PGFEM_fprintf(out,"</Points>\n");

  /* Elements */
  PGFEM_fprintf(out,"<Cells>\n");

  /* Connectivity */
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  long *nod;
  for (int i=0; i<ne; i++){
    nod = aloc1l (elem[i].toe);
    elemnodes (i,elem[i].toe,nod,elem);
    for (int j=0; j<elem[i].toe; j++){
      PGFEM_fprintf(out,"%10ld ",nod[j]);
    }
    PGFEM_fprintf(out,"\n");
    free(nod);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Nodes in element */
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  long off = 0;
  for (int i=0; i<ne; i++){
    off += elem[i].toe;
    PGFEM_fprintf(out,"%ld ",off);
  }
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Types */
  PGFEM_fprintf(out,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i=0; i<ne; i++){
    if(elem[i].toe == 4) PGFEM_fprintf(out,"10 ");
    else if(elem[i].toe == 8) PGFEM_fprintf(out,"12 ");
    else if(elem[i].toe == 10) PGFEM_fprintf(out,"24 ");
  }
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"</Cells>\n");


  /* Print Footer */
  PGFEM_fprintf(out,"</Piece>\n");
  PGFEM_fprintf(out,"</UnstructuredGrid>\n");
  PGFEM_fprintf(out,"</VTKFile>\n");

  /* Close file */
  fclose(out);
#endif
}


void VTK_print_cohesive_vtu(char *path,
          char *base_name,
          int time,
          int myrank,
          long nce,
          NODE *node,
          COEL *coel,
          SUPP sup,
          double *r,
          ENSIGHT ensight,
          const PGFem3D_opt *opts,
          const int mp_id)
{
  char cur_dir = '.';
  char *ptr_path = path;
  char dir_name[500];
  if(ptr_path == NULL) ptr_path = &cur_dir;
  sprintf(dir_name,"%s/%s/%s%.5d",ptr_path,out_dir,step_dir,time);
  if(make_path(dir_name,DIR_MODE) != 0){
    PGFEM_printf("Directory (%s) not created!\n",dir_name);
    abort();
  }

  /* Build filename and open file */
  char *filename;
  filename = (char*) PGFEM_calloc(500, sizeof(char));
  sprintf(filename,"%s/%s_cohesive_%d_%d.vtu",dir_name,base_name,myrank,time);
  FILE *out;
  if((out = fopen(filename,"w")) == NULL){
    PGFEM_printerr("ERROR: Cannot open file to write pvtu file! (%s)\n",
      filename);
    abort();
  }
  free(filename);

  /* Print header information */
  PGFEM_fprintf(out,"<?xml version=\"1.0\"?>\n");
  PGFEM_fprintf(out,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
    " byte_order=\"LittleEndian\">\n");
  PGFEM_fprintf(out,"<UnstructuredGrid>\n");
  PGFEM_fprintf(out,"<Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\" >\n",
    ensight->ncn,nce);

  /* POINT DATA */
  PGFEM_fprintf(out,"<PointData>\n");
  if(opts->multi_scale){
    PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"MacroDisplacement\""
      " NumberOfComponents=\"3\">\n");
    double u[3] = {0.0,0.0,0.0};
    for (int i=0; i<ensight->ncn;i++){
      /* total lagrangian, need only one of the nodes (coincide at time 0) */
      compute_interface_macro_disp_at_node(u,&node[ensight->Sm[i]],
             sup->F0,opts->analysis_type);
      PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",u[0],u[1],u[2]);
    }
    PGFEM_fprintf(out,"</DataArray>\n");
  }

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"Displacement\""
    " NumberOfComponents=\"3\">\n");
  for (int i=0;i<ensight->ncn;i++){
    for(int j=0; j<3; j++){
      if (node[ensight->Sm[i]].id_map[mp_id].id[j] == 0 && node[ensight->Sp[i]].id_map[mp_id].id[j] == 0)
  PGFEM_fprintf(out,"%12.12e ",0.0);

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] == 0 && node[ensight->Sp[i]].id_map[mp_id].id[j] >  0)
  PGFEM_fprintf(out,"%12.12e ",1./2.*(0.0 + r[node[ensight->Sp[i]].id_map[mp_id].id[j]-1]));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] == 0 && node[ensight->Sp[i]].id_map[mp_id].id[j] < 0)
  PGFEM_fprintf(out,"%12.12e ",
    1./2.*(0.0 + sup->defl[abs(node[ensight->Sp[i]].id_map[mp_id].id[j])-1]));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] >  0 && node[ensight->Sp[i]].id_map[mp_id].id[j] == 0)
  PGFEM_fprintf(out,"%12.12e ",1./2.*(r[node[ensight->Sm[i]].id_map[mp_id].id[j]-1] + 0.0));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] <  0 && node[ensight->Sp[i]].id_map[mp_id].id[j] == 0)
  PGFEM_fprintf(out,"%12.12e ",
    1./2.*(sup->defl[abs(node[ensight->Sm[i]].id_map[mp_id].id[j])-1] + 0.0));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] >  0 && node[ensight->Sp[i]].id_map[mp_id].id[j] <  0)
  PGFEM_fprintf(out,"%12.12e ",
    1./2.*(r[node[ensight->Sm[i]].id_map[mp_id].id[j]-1]
           + sup->defl[abs(node[ensight->Sp[i]].id_map[mp_id].id[j])-1]));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] <  0 && node[ensight->Sp[i]].id_map[mp_id].id[j] >  0)
  PGFEM_fprintf(out,"%12.12e ",
    1./2.*(sup->defl[abs(node[ensight->Sm[i]].id_map[mp_id].id[j])-1]
           + r[node[ensight->Sp[i]].id_map[mp_id].id[j]-1]));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] <  0 && node[ensight->Sp[i]].id_map[mp_id].id[j] <  0)
  PGFEM_fprintf(out,"%12.12e ",
    1./2.*(sup->defl[abs(node[ensight->Sm[i]].id_map[mp_id].id[j])-1]
           + sup->defl[abs(node[ensight->Sp[i]].id_map[mp_id].id[j])-1]));

      if (node[ensight->Sm[i]].id_map[mp_id].id[j] >  0 && node[ensight->Sp[i]].id_map[mp_id].id[j] >  0)
  PGFEM_fprintf(out,"%12.12e ",
    1./2.*(r[node[ensight->Sm[i]].id_map[mp_id].id[j]-1]
           + r[node[ensight->Sp[i]].id_map[mp_id].id[j]-1]));
    }
    PGFEM_fprintf(out,"\n");
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  PGFEM_fprintf(out,"</PointData>\n");

  /* CELL DATA */
  PGFEM_fprintf(out,"<CellData>\n");
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"Property\""
    " NumberOfComponents=\"1\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%ld\n",coel[i].pr);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"OpeningDisp\""
    " NumberOfComponents=\"3\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",
      coel[i].Xi[0],coel[i].Xi[1],coel[i].Xi[2]);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"Tractions\""
    " NumberOfComponents=\"3\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",
      coel[i].ti[0],coel[i].ti[1],coel[i].ti[2]);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"EffOpening\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e\n",coel[i].Xxi);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"NormalOpening\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e\n",coel[i].Xn);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"ShearOpening\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e\n",coel[i].Xs);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"EffTraction\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e\n",coel[i].txi);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"NormalTraction\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e\n",coel[i].tn);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Float64\" Name=\"ShearTraction\">\n");
  for(int i=0; i<nce; i++){
    PGFEM_fprintf(out,"%12.12e\n",coel[i].ts);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"</CellData>\n");

  /* Nodes */
  PGFEM_fprintf(out,"<Points>\n");
  PGFEM_fprintf(out,"<DataArray type=\"Float32\" NumberOfComponents=\"3\""
    " format=\"ascii\">\n");
  /* NOTE: tractions are in reference configuration, thus coordinates
     are in undeformed config. Can be warped using the Displacements
     vector */
  for(int i=0; i<ensight->ncn; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",
      0.5*(node[ensight->Sm[i]].x1_fd + node[ensight->Sp[i]].x1_fd),
      0.5*(node[ensight->Sm[i]].x2_fd + node[ensight->Sp[i]].x2_fd),
      0.5*(node[ensight->Sm[i]].x3_fd + node[ensight->Sp[i]].x3_fd));
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  PGFEM_fprintf(out,"</Points>\n");

  PGFEM_fprintf(out,"<Cells>\n");
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"connectivity\""
    " format=\"ascii\">\n");
  for (int i=0; i<nce; i++){
     for (int j=0; j<coel[i].toe/2; j++){
      PGFEM_fprintf(out,"%10ld ",ensight->No[coel[i].nod[j]]);
    }
    PGFEM_fprintf(out,"\n");
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"offsets\""
    " format=\"ascii\">\n");
  long off = 0;
  for (int i=0; i<nce; i++){
    off += coel[i].toe/2;
    PGFEM_fprintf(out,"%ld ",off);
  }
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i=0; i<nce; i++){
    if(coel[i].toe == 6) PGFEM_fprintf(out,"5 ");/* TRIANGLE */
    else if(coel[i].toe == 8) PGFEM_fprintf(out,"9 ");/* QUAD */
  }
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"</Cells>\n");

  /* Print Footer */
  PGFEM_fprintf(out,"</Piece>\n");
  PGFEM_fprintf(out,"</UnstructuredGrid>\n");
  PGFEM_fprintf(out,"</VTKFile>\n");

  /* Close file */
  fclose(out);
}

int VTK_get_filename(char *filename, 
                     const PGFem3D_opt *opts,
                     int time,
                     int myrank,
                     int is_it_master)
{
  int err = 0;
  
  char c_dir = '.';
  char *o_dir = "VTK";
  char *s_dir = "STEP_";  
  
  char *ptr_path = opts->opath;
  char *base_name = opts->ofname;
  
  char dir_name[1024];
  if(ptr_path == NULL) ptr_path = &c_dir;
  
  if(is_it_master)
    sprintf(dir_name,"%s/%s",ptr_path,o_dir);
  else
    sprintf(dir_name,"%s/%s/%s%.5d",ptr_path,o_dir,s_dir,time);

  if(make_path(dir_name,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",dir_name);
    abort();
  }
  
  // Build filename
  if(is_it_master)
    sprintf(filename,"%s/%s_%d.pvtu",dir_name,base_name,time);
  else
    sprintf(filename,"%s/%s_%d_%d.vtu",dir_name,base_name,myrank,time);

  return err;
}

int VTK_write_multiphysics_header(FILE *out, int is_it_master, long nodeno, long elemno)
{
  int err = 0;
  char vtk_file_format[1024];
  if(is_it_master)
    sprintf(vtk_file_format,"PUnstructuredGrid");
  else
    sprintf(vtk_file_format,"UnstructuredGrid");
    
  /* Print header information */
  PGFEM_fprintf(out,"<?xml version=\"1.0\"?>\n");
  PGFEM_fprintf(out,"<VTKFile type=\"%s\" version=\"0.1\""
                    " byte_order=\"LittleEndian\">\n", vtk_file_format);
  PGFEM_fprintf(out,"<%s>\n", vtk_file_format);

  if(is_it_master != 1)
    PGFEM_fprintf(out,"<Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\" >\n", nodeno, elemno);
  return err;
}

int VTK_write_multiphysics_master_footer(FILE *out,
                                         const PGFem3D_opt *opts,
                                         int nproc,
                                         int time)
{
  int err = 0;
  char *s_dir = "STEP_";
  char *base_name = opts->ofname;
  char dir_name[1024];  
  char filename[1024];
  /* write piece list */
  sprintf(dir_name,"%s%.5d",s_dir,time);
  for (int ia=0; ia< nproc; ia++){
    sprintf(filename,"%s/%s_%d_%d.vtu",dir_name,base_name,ia,time);
    PGFEM_fprintf(out,"<Piece Source=\"%s\"/>\n",filename);
  }
  /* write footer */
  PGFEM_fprintf(out,"</PUnstructuredGrid>\n");
  PGFEM_fprintf(out,"</VTKFile>\n");
  return err;
}

int VTK_write_multiphysics_DataArray_header(FILE *out,
                                            PRINT_MULTIPHYSICS_RESULT *pmr)
{
  int err = 0;

  PGFEM_fprintf(out,"<DataArray type=\"%s\" Name=\"%s\""
                    " NumberOfComponents=\"%d\" format=\"ascii\">\n", pmr->data_type,
                                                                      pmr->variable_name, 
                                                                      pmr->m_col);    
  return err;
}

int VTK_write_multiphysics_DataArray_footer(FILE *out)
{
  int err = 0;
  PGFEM_fprintf(out,"</DataArray>\n");
  return err;
}

int VTK_write_multiphysics_master(PRINT_MULTIPHYSICS_RESULT *pD,
                                  int datano,
                                  const PGFem3D_opt *opts,
                                  int time,
                                  int myrank,
                                  int nproc)
{
  int err = 0;
  char filename[1024];
  int is_it_master = 1;
  
  err += VTK_get_filename(filename,opts,time,myrank,is_it_master);
                     
  FILE *out = fopen(filename,"w");
  if(out == NULL)
  {
    PGFEM_printerr("ERROR: Cannot open file to write pvtu file! (%s)\n", filename);
    abort();
  }

  // write header information
  err += VTK_write_multiphysics_header(out,is_it_master,0,0);

  // write datatypes
  // Point Data
  PGFEM_fprintf(out,"<PPointData>\n");
  for(int ia=0; ia<datano; ia++)
  {
    if(pD[ia].is_point_data) // point data
    {
      PGFEM_fprintf(out,"<PDataArray type=\"%s\" Name=\"%s\""
                        " NumberOfComponents=\"%d\"/>\n", pD[ia].data_type,
                                                          pD[ia].variable_name,
                                                          pD[ia].m_col);
    }
  }
  PGFEM_fprintf(out,"</PPointData>\n");
  
  // Cell data
  PGFEM_fprintf(out,"<PCellData>\n");
  for(int ia=0; ia<datano; ia++)
  {
    if(pD[ia].is_point_data != 1) // cell data
    {
      PGFEM_fprintf(out,"<PDataArray type=\"%s\" Name=\"%s\""
                        " NumberOfComponents=\"%d\"/>\n", pD[ia].data_type,
                                                          pD[ia].variable_name,
                                                          pD[ia].m_col);
    }
  }
  PGFEM_fprintf(out,"</PCellData>\n");  

  /* Points */
  PGFEM_fprintf(out,"<PPoints>\n");
  PGFEM_fprintf(out,"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
  PGFEM_fprintf(out,"</PPoints>\n");


  err += VTK_write_multiphysics_master_footer(out,opts,nproc,time);

  // close file
  fclose(out);
  return err;
}

int VTK_write_mesh(GRID *grid,
                   FILE *out)
{
  int err = 0;
  /* Nodes */
  NODE    *node = grid->node;
  ELEMENT *elem = grid->element;
  
  PGFEM_fprintf(out,"<Points>\n");
  PGFEM_fprintf(out,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0; i<grid->nn; i++){
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",
      node[i].x1_fd,node[i].x2_fd,node[i].x3_fd);
  }
  PGFEM_fprintf(out,"</DataArray>\n");
  PGFEM_fprintf(out,"</Points>\n");

  /* Elements */
  PGFEM_fprintf(out,"<Cells>\n");

  /* Connectivity */
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  long *nod;
  for (int i=0; i<grid->ne; i++){
    nod = aloc1l (elem[i].toe);
    elemnodes (i,elem[i].toe,nod,elem);
    for (int j=0; j<elem[i].toe; j++){
      PGFEM_fprintf(out,"%10ld ",nod[j]);
    }
    PGFEM_fprintf(out,"\n");
    free(nod);
  }
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Nodes in element */
  PGFEM_fprintf(out,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  long off = 0;
  for (int i=0; i<grid->ne; i++){
    off += elem[i].toe;
    PGFEM_fprintf(out,"%ld ",off);
  }
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"</DataArray>\n");

  /* Types */
  PGFEM_fprintf(out,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i=0; i<grid->ne; i++){
    if(elem[i].toe == 4) PGFEM_fprintf(out,"10 ");
    else if(elem[i].toe == 8) PGFEM_fprintf(out,"12 ");
    else if(elem[i].toe == 10) PGFEM_fprintf(out,"24 ");
  }
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"</DataArray>\n");

  PGFEM_fprintf(out,"</Cells>\n");
  
  return err;  
}


int VTK_write_multiphysics_vtu(GRID *grid,
                               FIELD_VARIABLES *FV,
                               LOADING_STEPS *load,
                               PRINT_MULTIPHYSICS_RESULT *pD,
                               int datano,
                               const PGFem3D_opt *opts,
                               int time,
                               int myrank)
{
  int err = 0;
  char filename[1024];
  int is_it_master = 0;
  err += VTK_get_filename(filename, opts, time, myrank, 0);

  FILE *out = fopen(filename,"w");
  if(out == NULL)
  {
    PGFEM_printerr("ERROR: Cannot open file to write pvtu file! (%s)\n", filename);
    abort();
  }
  
  err += VTK_write_multiphysics_header(out,is_it_master,grid->nn,grid->ne);

  PGFEM_fprintf(out,"<PointData>\n");   
  for(int ia=0; ia<datano; ia++)
  {  
    if(pD[ia].is_point_data) // point data
      err += pD[ia].write_vtk(out,grid,FV,load,pD+ia,opts);
  }
  PGFEM_fprintf(out,"</PointData>\n"); 
  
  PGFEM_fprintf(out,"<CellData>\n"); 
  for(int ia=0; ia<datano; ia++)
  {
    if(pD[ia].is_point_data != 1) // Cell data
      err += pD[ia].write_vtk(out,grid,FV,load,pD+ia,opts); 
  }
  PGFEM_fprintf(out,"</CellData>\n");   
  
  err += VTK_write_mesh(grid, out);
  
  PGFEM_fprintf(out,"</Piece>\n");  
  PGFEM_fprintf(out,"</UnstructuredGrid>\n");
  PGFEM_fprintf(out,"</VTKFile>\n");
  
  fclose(out);
  return err;
}

int VTK_write_data_double(FILE *out,
                          GRID *grid,
                          FIELD_VARIABLES *FV,                          
                          LOADING_STEPS *load,
                          PRINT_MULTIPHYSICS_RESULT *pmr,
                          const PGFem3D_opt *opts)
{
  int err = 0;
  double *p_data = (double *) pmr->p_data;
  
  err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int ia=0; ia<pmr->m_row; ia++)
  {
    for (int ja=0; ja<pmr->m_col; ja++)
      PGFEM_fprintf(out,"%12.12e ", p_data[ia*(pmr->m_col)+ja]);

    PGFEM_fprintf(out,"\n");
  }
  err += VTK_write_multiphysics_DataArray_footer(out);
  return err;
}

int VTK_write_data_MacroDisplacement(FILE *out,
                                     GRID *grid,
                                     FIELD_VARIABLES *FV,                          
                                     LOADING_STEPS *load,
                                     PRINT_MULTIPHYSICS_RESULT *pmr,
                                     const PGFem3D_opt *opts)
{
  int err = 0;
  
  NODE *node = grid->node;
  SUPP sup = load->sups[pmr->mp_id];  

  if(sup->multi_scale){
      err += VTK_write_multiphysics_DataArray_header(out, pmr);  
    double *jump = aloc1(3);
    compute_macro_grad_u(sup->F0,sup,opts->analysis_type);
    for(int i=0; i<grid->nn; i++){
      compute_interface_macro_disp_at_node(jump,&node[i],sup->F0,opts->analysis_type);
      PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",jump[0],jump[1],jump[2]);
    }
    err += VTK_write_multiphysics_DataArray_footer(out);
    free(jump);
  }
  
  return err;
}

int VTK_write_data_Nodal_pressure(FILE *out,
                                  GRID *grid,
                                  FIELD_VARIABLES *FV,                          
                                  LOADING_STEPS *load,
                                  PRINT_MULTIPHYSICS_RESULT *pmr,
                                  const PGFem3D_opt *opts)
{
  int err = 0;
  ELEMENT *elem = grid->element;
  NODE *node = grid->node;
  int mp_id = pmr->mp_id;
  /* Pressure */
  switch(opts->analysis_type){
  case TF:
    if(elem[0].toe==10)
      break;
  case STABILIZED:
  case MINI:
  case MINI_3F:
      err += VTK_write_multiphysics_DataArray_header(out, pmr);  
    for (int i=0; i<grid->nn; i++){
      int id = node[i].id_map[mp_id].id[3]-1;
      PGFEM_fprintf(out,"%12.12e\n",FV[mp_id].u_np1[id]);
    }
    err += VTK_write_multiphysics_DataArray_footer(out);
    break;
  default:
    break;
  }
  return err;
}

int VTK_write_data_CauchyStress(FILE *out,
                                GRID *grid,
                                FIELD_VARIABLES *FV,                          
                                LOADING_STEPS *load,
                                PRINT_MULTIPHYSICS_RESULT *pmr,
                                const PGFem3D_opt *opts)
{
  int err = 0;
  
  SIG *sig = FV[pmr->mp_id].sig;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int ia=0; ia<grid->ne; ia++)
  {
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e\n",
                  sig[ia].el.o[0],sig[ia].el.o[1],sig[ia].el.o[2],
                  sig[ia].el.o[5],sig[ia].el.o[3],sig[ia].el.o[4]);
  }
  err += VTK_write_multiphysics_DataArray_footer(out);    
  return err;
}

int VTK_write_data_EulerStrain(FILE *out,
                               GRID *grid,
                               FIELD_VARIABLES *FV,                          
                               LOADING_STEPS *load,
                               PRINT_MULTIPHYSICS_RESULT *pmr,
                               const PGFem3D_opt *opts)
{
  int err = 0;
  
  EPS *eps = FV[pmr->mp_id].eps;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
  {
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e\n",
                  eps[i].el.o[0],eps[i].el.o[1],eps[i].el.o[2],
                  eps[i].el.o[5]/2,eps[i].el.o[3]/2,eps[i].el.o[4]/2);
  }
  err += VTK_write_multiphysics_DataArray_footer(out); 
  return err;
}  

int VTK_write_data_EffectiveStrain(FILE *out,
                                   GRID *grid,
                                   FIELD_VARIABLES *FV,                          
                                   LOADING_STEPS *load,
                                   PRINT_MULTIPHYSICS_RESULT *pmr,
                                   const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
    PGFEM_fprintf(out,"%12.12e\n",eps[i].el.eq);
  
  err += VTK_write_multiphysics_DataArray_footer(out); 
  return err;
}

int VTK_write_data_EffectiveStress(FILE *out,
                               GRID *grid,
                               FIELD_VARIABLES *FV,                          
                               LOADING_STEPS *load,
                               PRINT_MULTIPHYSICS_RESULT *pmr,
                               const PGFem3D_opt *opts)
{
  int err = 0;
  SIG *sig = FV[pmr->mp_id].sig;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
    PGFEM_fprintf(out,"%12.12e\n",sig[i].el.eq);
  
  err += VTK_write_multiphysics_DataArray_footer(out); 
  return err;
}

int VTK_write_data_CellProperty(FILE *out,
                                GRID *grid,
                                FIELD_VARIABLES *FV,                          
                                LOADING_STEPS *load,
                                PRINT_MULTIPHYSICS_RESULT *pmr,
                                const PGFem3D_opt *opts)
{
  int err = 0;
  ELEMENT *elem = grid->element;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
    PGFEM_fprintf(out,"%ld\n",elem[i].pr);
  
  return err;
}


int VTK_write_data_Damage(FILE *out,
                          GRID *grid,
                          FIELD_VARIABLES *FV,                          
                          LOADING_STEPS *load,
                          PRINT_MULTIPHYSICS_RESULT *pmr,
                          const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
    PGFEM_fprintf(out,"%12.12e\n",eps[i].dam[0].wn); // supported for Linear element (1 ip)
  
  err += VTK_write_multiphysics_DataArray_footer(out); 
  return err;
}

int VTK_write_data_Chi(FILE *out,
                       GRID *grid,
                       FIELD_VARIABLES *FV,                          
                       LOADING_STEPS *load,
                       PRINT_MULTIPHYSICS_RESULT *pmr,
                       const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
    PGFEM_fprintf(out,"%12.12e\n",eps[i].dam[0].Xn); // supported for Linear element (1 ip)
  
  err += VTK_write_multiphysics_DataArray_footer(out); 
  return err;
}

int VTK_write_data_F(FILE *out,
                     GRID *grid,
                     FIELD_VARIABLES *FV,                          
                     LOADING_STEPS *load,
                     PRINT_MULTIPHYSICS_RESULT *pmr,
                     const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
  {
    const double *F = eps[i].il[0].F;
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e "
	        	          "%12.12e %12.12e %12.12e %12.12e\n",
		                  F[0],F[1],F[2],F[3],F[4],F[5],F[6],F[7],F[8]);
  }    
  
  return err;
}

int VTK_write_data_P(FILE *out,
                     GRID *grid,
                     FIELD_VARIABLES *FV,                          
                     LOADING_STEPS *load,
                     PRINT_MULTIPHYSICS_RESULT *pmr,
                     const PGFem3D_opt *opts)
{
  int err = 0;
  SIG *sig = FV[pmr->mp_id].sig;
  EPS *eps = FV[pmr->mp_id].eps;
  double S[9], P[9];
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int i=0; i<grid->ne; i++)
  {
    const double *F = eps[i].il[0].F;
    S[0] = sig[i].il[0].o[0];
    S[1] = sig[i].il[0].o[5];
    S[2] = sig[i].il[0].o[4];
    
    S[3] = sig[i].il[0].o[5];
    S[4] = sig[i].il[0].o[1];
    S[5] = sig[i].il[0].o[3];

    S[6] = sig[i].il[0].o[3];
    S[7] = sig[i].il[0].o[4];
    S[8] = sig[i].il[0].o[2];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ndim,ndim,ndim,
                1.0,F,ndim,S,ndim,0.0,P,ndim);    
    
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e %12.12e "
	        	          "%12.12e %12.12e %12.12e %12.12e\n",
		                  P[0],P[1],P[2],P[3],P[4],P[5],P[6],P[7],P[8]);
  }    
  err += VTK_write_multiphysics_DataArray_footer(out);
  return err;
}

int VTK_write_data_W(FILE *out,
                     GRID *grid,
                     FIELD_VARIABLES *FV,                          
                     LOADING_STEPS *load,
                     PRINT_MULTIPHYSICS_RESULT *pmr,
                     const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  double S[9], P[9];
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for(int i=0; i<grid->ne; i++)
    PGFEM_fprintf(out,"%12.12e\n",eps[i].il[0].Y);
  
  err += VTK_write_multiphysics_DataArray_footer(out); 
  return err;
}    

int VTK_write_data_ElementPressure(FILE *out,
                                  GRID *grid,
                                  FIELD_VARIABLES *FV,                          
                                  LOADING_STEPS *load,
                                  PRINT_MULTIPHYSICS_RESULT *pmr,
                                  const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  ELEMENT *elem = grid->element;
  double S[9], P[9];
  
  if(opts->analysis_type==TF && elem[0].toe==10)
  {  
      err += VTK_write_multiphysics_DataArray_header(out, pmr);  
    for (int i=0; i<grid->ne; i++)
      PGFEM_fprintf(out,"%12.12e\n",eps[i].d_T[0]);
  
    err += VTK_write_multiphysics_DataArray_footer(out); 
  }
  return err;
}

int VTK_write_data_ElementVolume(FILE *out,
                                  GRID *grid,
                                  FIELD_VARIABLES *FV,                          
                                  LOADING_STEPS *load,
                                  PRINT_MULTIPHYSICS_RESULT *pmr,
                                  const PGFem3D_opt *opts)
{
  int err = 0;
  EPS *eps = FV[pmr->mp_id].eps;
  double S[9], P[9];
  
  if(opts->analysis_type==TF)
  {  
      err += VTK_write_multiphysics_DataArray_header(out, pmr);  
    for (int i=0; i<grid->ne; i++)
      PGFEM_fprintf(out,"%12.12e\n",eps[i].T[0]);
  
    err += VTK_write_multiphysics_DataArray_footer(out); 
  }
  return err;
}

int VTK_write_data_HeatFlux(FILE *out,
                            GRID *grid,
                            FIELD_VARIABLES *FV,                          
                            LOADING_STEPS *load,
                            PRINT_MULTIPHYSICS_RESULT *pmr,
                            const PGFem3D_opt *opts)
{
  int err = 0;
  
  SIG *sig = FV[pmr->mp_id].sig;
  
    err += VTK_write_multiphysics_DataArray_header(out, pmr);  
  for (int ia=0; ia<grid->ne; ia++)
  {
    PGFEM_fprintf(out,"%12.12e %12.12e %12.12e\n",
                  sig[ia].el.o[0],sig[ia].el.o[1],sig[ia].el.o[2]);
  }
  err += VTK_write_multiphysics_DataArray_footer(out);    
  return err;
}

int VTK_construct_PMR(GRID *grid,
                      FIELD_VARIABLES *FV,
                      MULTIPHYSICS *mp,
                      PRINT_MULTIPHYSICS_RESULT *pmr)
{
  int cnt_pmr = 0;
  for(int ia=0; ia<mp->physicsno; ia++)
  {
    int physics_id = mp->physics_ids[ia];
    switch(physics_id)
    {
      case MULTIPHYSICS_MECHANICAL:
      {
        for(int ib=0; ib<mp->write_no[ia]; ib++)
        {
          pmr[cnt_pmr].is_point_data = 0;          
          pmr[cnt_pmr].mp_id      = ia;
          pmr[cnt_pmr].physics_id = MULTIPHYSICS_MECHANICAL;
          pmr[cnt_pmr].m_row      = grid->ne;
          pmr[cnt_pmr].m_col      = 1;
          pmr[cnt_pmr].p_data     = NULL; 
          sprintf(pmr[cnt_pmr].data_type, "Float64");
          switch(mp->write_ids[ia][ib])
          {
            case MECHANICAL_Var_Displacement:
              pmr[cnt_pmr].is_point_data = 1;
              pmr[cnt_pmr].m_row         = grid->nn;
              pmr[cnt_pmr].m_col         = 3;
              pmr[cnt_pmr].p_data        = FV[ia].u_n;        
              pmr[cnt_pmr].write_vtk     = VTK_write_data_double;
              sprintf(pmr[cnt_pmr].variable_name, "Displacement");
              break;
            case MECHANICAL_Var_MacroDisplacement:
              pmr[cnt_pmr].is_point_data = 1;
              pmr[cnt_pmr].m_row         = grid->nn;
              pmr[cnt_pmr].m_col         = 3;
              pmr[cnt_pmr].write_vtk = VTK_write_data_MacroDisplacement;
              sprintf(pmr[cnt_pmr].variable_name, "MacroDisplacement");
              break;
            case MECHANICAL_Var_NodalPressure:
              pmr[cnt_pmr].is_point_data = 1;
              pmr[cnt_pmr].m_row         = grid->nn;
              pmr[cnt_pmr].m_col         = 1;
              pmr[cnt_pmr].write_vtk     = VTK_write_data_Nodal_pressure;
              sprintf(pmr[cnt_pmr].variable_name, "Pressure");
              break;
            case MECHANICAL_Var_CauchyStress:
              pmr[cnt_pmr].m_col     = 6;
              pmr[cnt_pmr].write_vtk = VTK_write_data_CauchyStress;
              sprintf(pmr[cnt_pmr].variable_name, "CauchyStress");
              break;
            case MECHANICAL_Var_EulerStrain:
              pmr[cnt_pmr].m_col     = 6;
              pmr[cnt_pmr].write_vtk = VTK_write_data_EulerStrain;
              sprintf(pmr[cnt_pmr].variable_name, "EulerStrain");
              break;
            case MECHANICAL_Var_EffectiveStrain:
              pmr[cnt_pmr].write_vtk = VTK_write_data_EffectiveStrain;
              sprintf(pmr[cnt_pmr].variable_name, "EffectiveStrain");
              break;
            case MECHANICAL_Var_EffectiveStress:
              pmr[cnt_pmr].write_vtk = VTK_write_data_EffectiveStress;
              sprintf(pmr[cnt_pmr].variable_name, "EffectiveStress");
              break;
            case MECHANICAL_Var_CellProperty:
              pmr[cnt_pmr].write_vtk = VTK_write_data_CellProperty;
              sprintf(pmr[cnt_pmr].data_type, "Int64");
              sprintf(pmr[cnt_pmr].variable_name, "CellProperty");
              break;
            case MECHANICAL_Var_Damage:
              pmr[cnt_pmr].write_vtk = VTK_write_data_Damage;
              sprintf(pmr[cnt_pmr].variable_name, "Damage");
              break;
            case MECHANICAL_Var_Chi:
              pmr[cnt_pmr].write_vtk = VTK_write_data_Chi;
              sprintf(pmr[cnt_pmr].variable_name, "Chi");
              break;
            case MECHANICAL_Var_F:
              pmr[cnt_pmr].m_col     = 9;
              pmr[cnt_pmr].write_vtk = VTK_write_data_F;
              sprintf(pmr[cnt_pmr].variable_name, "F");
              break;
            case MECHANICAL_Var_P:
              pmr[cnt_pmr].m_col     = 9;
              pmr[cnt_pmr].write_vtk = VTK_write_data_P;
              sprintf(pmr[cnt_pmr].variable_name, "P");
              break;
            case MECHANICAL_Var_W:
              pmr[cnt_pmr].write_vtk = VTK_write_data_W;
              sprintf(pmr[cnt_pmr].variable_name, "W");
              break;
            case MECHANICAL_Var_ElementPressure:
              pmr[cnt_pmr].write_vtk = VTK_write_data_ElementPressure;
              sprintf(pmr[cnt_pmr].variable_name, "TF_Pressure");
              break;
            case MECHANICAL_Var_ElementVolume:
              pmr[cnt_pmr].write_vtk = VTK_write_data_ElementVolume;
              sprintf(pmr[cnt_pmr].variable_name, "TF_Volume");
              break;
            default:
              pmr[cnt_pmr].is_point_data = 1;
              pmr[cnt_pmr].m_row         = grid->nn;
              pmr[cnt_pmr].m_col         = 3;
              pmr[cnt_pmr].p_data        = FV[ia].u_n;        
              pmr[cnt_pmr].write_vtk     = VTK_write_data_double;
              sprintf(pmr[cnt_pmr].variable_name, "Displacement");              
              break;
          }
          cnt_pmr++;
        }
        break;
      }
      case MULTIPHYSICS_THERMAL:
      {
        for(int ib=0; ib<mp->write_no[ia]; ib++)
        {
          pmr[cnt_pmr].is_point_data = 0;          
          pmr[cnt_pmr].mp_id      = ia;
          pmr[cnt_pmr].physics_id = MULTIPHYSICS_THERMAL;
          pmr[cnt_pmr].m_row      = grid->ne;
          pmr[cnt_pmr].m_col      = 1;
          pmr[cnt_pmr].p_data     = NULL; 
          sprintf(pmr[cnt_pmr].data_type, "Float64");

          switch(mp->write_ids[ia][ib])
          {
            case THERMAL_Var_Temperature:
              pmr[cnt_pmr].is_point_data = 1;
              pmr[cnt_pmr].m_row         = grid->nn;
              pmr[cnt_pmr].m_col         = 1;
              pmr[cnt_pmr].p_data        = FV[ia].u_n;        
              pmr[cnt_pmr].write_vtk     = VTK_write_data_double;
              sprintf(pmr[cnt_pmr].variable_name, "Temperature");
              break;
            case THERMAL_Var_HeatFlux:
              pmr[cnt_pmr].m_col         = 3;
              pmr[cnt_pmr].write_vtk     = VTK_write_data_HeatFlux;
              sprintf(pmr[cnt_pmr].variable_name, "HeatFlux");
              break;
            default:
              pmr[cnt_pmr].is_point_data = 1;
              pmr[cnt_pmr].m_row         = grid->nn;
              pmr[cnt_pmr].m_col         = 1;
              pmr[cnt_pmr].p_data        = FV[ia].u_n;        
              pmr[cnt_pmr].write_vtk     = VTK_write_data_double;
              sprintf(pmr[cnt_pmr].variable_name, "Temperature");
              break;                            
          }
          cnt_pmr++;           
        }                  
        break;
      }  
      case MULTIPHYSICS_CHEMICAL:
      {
        break;
      }
      default:
      {
      }
    }
  }  

  return 0;
}                                                   