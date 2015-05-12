/* HEADER */
#include "in.h"

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

SUPP read_supports (FILE *in,
		    long nn,
		    long ndofn,
		    NODE *node)
/*
  in    - Input file
  nl    - Number of layers
  ndofn - Number of degrees of freedom in one node
  node  - Structure type of NODE
  sup   - Structure type of SUPP
  
  %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
*/
{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if (PFEM_DEBUG) PGFEM_printf("[%d] reading supports.\n",err_rank);
  long i,k,ii,n,pom;
  SUPP sup;
  
  /************************/
  /* read supported nodes */
  /************************/
  
  sup = (SUPP) PGFEM_calloc (1,sizeof(SUPP_1));
  
  fscanf (in,"%ld",&sup->nsn);
  
  if (sup->nsn == 0)  sup->supp = (long *) PGFEM_calloc (1,sizeof(long));
  else                sup->supp = (long *) PGFEM_calloc (sup->nsn,sizeof(long));
  
  for (i=0;i<sup->nsn;i++){
    fscanf (in,"%ld",&n);
    sup->supp[i] = n;
    
    pom = 0;
    for (k=0;k<ndofn;k++){
      fscanf (in,"%ld",&node[n].id[k]);
      if (node[n].id[k] <= -1 && pom == 0) {
	sup->ndn++; pom = 1;
      }
    }/* end k */
    if(ferror(in)){
      PGFEM_printerr("[%d]ERROR:fscanf returned error"
	      " reading support %ld!\n",err_rank,i);
      PGFEM_Abort();
    } else if(feof(in)){
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
	      err_rank);
      PGFEM_Abort();
    }
  }/* end i */
  
  /***************************************************/
  /* create list of nodes with prescribed deflection */
  /***************************************************/
  
  if (sup->ndn == 0)  sup->lnpd = (long *) PGFEM_calloc (1,sizeof(long));
  else                sup->lnpd = (long *) PGFEM_calloc (sup->ndn,sizeof(long));
  
  ii = 0;
  for (i=0;i<sup->nsn;i++){
    for (k=0;k<ndofn;k++){
      if (node[sup->supp[i]].id[k] <= -1) {
	sup->lnpd[ii] = sup->supp[i];
	ii++;  break;
      }
    }/* end k */
  }/* end i */
  
  /*****************************************/
  /* read nodes with prescribed deflection */
  /*****************************************/
  
  fscanf (in,"%ld",&sup->npd);
  
  if (sup->npd == 0){
    sup->defl   = (double *) PGFEM_calloc (1,sizeof(double));
    sup->defl_d = (double *) PGFEM_calloc (1,sizeof(double));
  }
  else{
    sup->defl   = (double *) PGFEM_calloc (sup->npd,sizeof(double));
    sup->defl_d = (double *) PGFEM_calloc (sup->npd,sizeof(double));
  }
  
  for (i=0;i<sup->npd;i++) fscanf (in,"%lf",&sup->defl_d[i]);

  if(ferror(in)){
    PGFEM_printerr("[%d]ERROR:fscanf returned error"
	    " reading prescribed deflections!\n",err_rank);
    PGFEM_Abort();
  } else if(feof(in)){
    PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
	    err_rank);
    PGFEM_Abort();
  }

  /* allocate the prescribed macro deformation gradient */
  sup->F0 = (double *) PGFEM_calloc(9,sizeof(double));
  sup->N0 = (double *) PGFEM_calloc(3,sizeof(double));

  return (sup);
}

int read_material (FILE *in,
                   const size_t mat_id,
                   MATERIAL *mater,
                   const int legacy)
{
  int err = 0;
  if (legacy) {
    fscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &mater[mat_id].Ex,&mater[mat_id].Ey,&mater[mat_id].Ez,
            &mater[mat_id].Gyz,&mater[mat_id].Gxz,&mater[mat_id].Gxy,
            &mater[mat_id].nyz,&mater[mat_id].nxz,&mater[mat_id].nxy,
            &mater[mat_id].ax,&mater[mat_id].ay,&mater[mat_id].az,
            &mater[mat_id].sig);
    mater[mat_id].devPotFlag = mater[mat_id].volPotFlag = 0;
  } else {
    fscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
            &mater[mat_id].Ex,&mater[mat_id].Ey,&mater[mat_id].Ez,
            &mater[mat_id].Gyz,&mater[mat_id].Gxz,&mater[mat_id].Gxy,
            &mater[mat_id].nyz,&mater[mat_id].nxz,&mater[mat_id].nxy,
            &mater[mat_id].ax,&mater[mat_id].ay,&mater[mat_id].az,
            &mater[mat_id].sig,&mater[mat_id].devPotFlag,&mater[mat_id].volPotFlag);
  }

  if(ferror(in)){
    PGFEM_printerr("ERROR: fscanf returned error in: %s(%s)\n",__func__,__FILE__);
    err++;
  } else if(feof(in)){
    PGFEM_printerr("ERROR: prematurely reached end of file in %s(%s)\n",__func__,__FILE__);
    err++;
  }
  return err;
}

void read_matgeom (FILE *in,
		   long nc,
		   long np,
		   MATGEOM matgeom)
/*

*/
{
  if (PFEM_DEBUG) PGFEM_printf("[%d] reading material geom.\n",1);
  long i,j;
  
  for (i=0;i<nc;i++){
    fscanf (in,"%lf",&matgeom->cf[i]);
    matgeom->cm[i] = 1.0 - matgeom->cf[i];
    matgeom->cd[i] = 0.0;
  }
  for (i=0;i<np;i++){
    for (j=0;j<9;j++){
      fscanf (in,"%lf",&matgeom->ee[i][j]);
    }
  }
  fscanf (in,"%ld %lf %lf",&matgeom->SH,&matgeom->a1,&matgeom->a2);

}

void read_nodal_load (FILE *in,
		      long nln,
		      long ndofn,
		      ZATNODE *znod)
/*
  
*/
{
  if (PFEM_DEBUG) PGFEM_printf("[%d] reading nodal loads.\n",1);
  long i,j;
  
  for (i=0;i<nln;i++){
    
    fscanf (in,"%ld",&znod[i].nod);
    
    for (j=0;j<ndofn;j++)
      fscanf (in,"%lf",&znod[i].load[j]);
  } 
 
}

void read_elem_surface_load (FILE *in,
			     long nle_s,
			     long ndofn,
			     ELEMENT *elem,
			     ZATELEM *zele_s)
     /*
       
     */
{
  long i,j,ii,jj;
  
  for (i=0;i<nle_s;i++){
    fscanf (in,"%ld",&zele_s[i].elem);
    ii = elem[zele_s[i].elem].toe;
    if (ii == 4)  jj = 3;
    if (ii == 8)  jj = 4;
    if (ii == 10) jj = 6;
    zele_s[i].sur = (long*) PGFEM_calloc (jj,sizeof(long));
    for (j=0;j<jj;j++){
      fscanf (in,"%ld",&zele_s[i].sur[j]);
    }
    for (j=0;j<ndofn;j++){
      fscanf (in,"%lf",&zele_s[i].load[j]);
    }
  }

}

int override_prescribed_displacements(SUPP sup,
				      const PGFem3D_opt *opt)
{
  int err = 0;
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if(err_rank == 0){
    PGFEM_printf("Overriding the prescribed displacements with:\n"
	   "%s\n",opt->pre_disp_file);
  }

  FILE *in = PGFEM_fopen(opt->pre_disp_file,"r");
  for(int i=0; i<sup->npd; i++){
    fscanf(in,"%lf",&sup->defl_d[i]);
  }
  err = ferror(in);
  PGFEM_fclose(in);
  return err;
}
