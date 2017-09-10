#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFEM_mpi.h"
#include "allocation.h"
#include "in.h"
#include "utils.h"
#include <cstring>
#include <cassert>

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

/// read Dirichlet boundary conditions on nodes
///
/// This function will generate and return SUPP object.
/// if in==NULL, assumed no boundary conditions. Code will proceed to next
/// step in reading boundary condition values.
///
/// \param[in]  in    Input file
/// \param[in]  ndofn Number of degrees of freedom in one node
/// \param[in]  node  Structure type of NODE
/// \param[in]  mp_id multiphysics id
/// \return     sup   created boundary condition structure
SUPP read_Dirichlet_BCs(FILE *in,
                        long nn,
                        long ndofn,
                        NODE *node,
                        const int mp_id)
{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if (PFEM_DEBUG) PGFEM_printf("[%d] reading supports.\n",err_rank);
  long n,pom;
  SUPP sup = PGFEM_calloc (SUPP_1, 1);

  // read supported nodes

  if(in==NULL)
    sup->nsn = 0;
  else
    CHECK_SCANF(in, "%ld", &sup->nsn);

  if(sup->nsn == 0)  sup->supp = PGFEM_calloc (long, 1);
  else               sup->supp = PGFEM_calloc (long, sup->nsn);

  for(int i=0; i<sup->nsn; i++)
  {
    CHECK_SCANF(in,"%ld",&n);
    sup->supp[i] = n;

    pom = 0;
    for(int k=0; k<ndofn; k++)
    {
      CHECK_SCANF(in,"%ld",&node[n].id_map[mp_id].id[k]);
      if((node[n].id_map[mp_id].id[k] == 1 || node[n].id_map[mp_id].id[k] <= -1) && pom == 0) {
        sup->ndn++; pom = 1;
      }
    }
    if(ferror(in))
    {
      PGFEM_printerr("[%d]ERROR:fscanf returned error"
                     " reading support %ld!\n",err_rank,i);
      PGFEM_Abort();
    }
    else if(feof(in))
    {
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                     err_rank);
      PGFEM_Abort();
    }
  }

  /***************************************************/
  /* create list of nodes with prescribed deflection */
  /***************************************************/

  if (sup->ndn == 0)  sup->lnpd = PGFEM_calloc (long, 1);
  else                sup->lnpd = PGFEM_calloc (long, sup->ndn);

  int ii = 0;
  for(int i=0; i<sup->nsn; i++)
  {
    for(int k=0; k<ndofn; k++)
    {
      if(node[n].id_map[mp_id].id[k] == 1 || node[n].id_map[mp_id].id[k] <= -1) {
        sup->lnpd[ii] = sup->supp[i];
        ii++;  break;
      }
    }
  }

  // allocate the prescribed macro deformation gradient
  sup->F0 = PGFEM_calloc(double, 9);
  sup->N0 = PGFEM_calloc(double, 3);

  return (sup);
}

/// read Dirichlet boundary condition values
///
/// Before read BC values, SUPP object should be created.
/// Prior to call this function, read_Dirichlet_BCs function should be called first.
///
/// \param[in]  in    Input file
/// \param[in]  ndofn Number of degrees of freedom in one node
/// \param[in]  node  Structure type of NODE
/// \param[in]  mp_id multiphysics id
/// \return non-zero on internal ERROR
int read_Dirichlet_BCs_values(FILE *in,
                              long nn,
                              long ndofn,
                              NODE *node,
                              SUPP sup,
                              const int mp_id)
{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if (PFEM_DEBUG) PGFEM_printf("[%d] reading BCs_values.\n",err_rank);

  // read nodes with prescribed deflection
  CHECK_SCANF(in,"%ld",&sup->npd);

  if (sup->npd == 0)
  {
    sup->defl   = PGFEM_calloc (double, 1);
    sup->defl_d = PGFEM_calloc (double, 1);
  }
  else
  {
    sup->defl   = PGFEM_calloc (double, sup->npd);
    sup->defl_d = PGFEM_calloc (double, sup->npd);
  }

  for (int i=0; i<sup->npd; i++) {
    CHECK_SCANF(in,"%lf",&sup->defl_d[i]);
  }

  if (ferror(in)) {
    PGFEM_printerr("[%d]ERROR:fscanf returned error"
                   " reading prescribed deflections!\n",err_rank);
    PGFEM_Abort();
  } else if(feof(in)){
    PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                   err_rank);
    PGFEM_Abort();
  }

  return 0;
}

SUPP read_supports(FILE *in,
                   long nn,
                   long ndofn,
                   NODE *node,
                   const int mp_id)
/*
  in    - Input file
  nl    - Number of layers
  ndofn - Number of degrees of freedom in one node
  node  - Structure type of NODE
  sup   - Structure type of SUPP

  %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
*/
{
  SUPP sup = read_Dirichlet_BCs(in,nn,ndofn,node,mp_id);
  read_Dirichlet_BCs_values(in,nn,ndofn,node,sup,mp_id);

  return sup;
}



int read_material (FILE *in,
                   const size_t mat_id,
                   Material *mater,
                   const int legacy)
{
  int err = 0;
  if (legacy) {
    CHECK_SCANF(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &mater[mat_id].Ex,&mater[mat_id].Ey,&mater[mat_id].Ez,
                &mater[mat_id].Gyz,&mater[mat_id].Gxz,&mater[mat_id].Gxy,
                &mater[mat_id].nyz,&mater[mat_id].nxz,&mater[mat_id].nxy,
                &mater[mat_id].ax,&mater[mat_id].ay,&mater[mat_id].az,
                &mater[mat_id].sig);

    mater[mat_id].devPotFlag = mater[mat_id].volPotFlag = 0;
  } else {
    CHECK_SCANF(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
                &mater[mat_id].Ex,&mater[mat_id].Ey,&mater[mat_id].Ez,
                &mater[mat_id].Gyz,&mater[mat_id].Gxz,&mater[mat_id].Gxy,
                &mater[mat_id].nyz,&mater[mat_id].nxz,&mater[mat_id].nxy,
                &mater[mat_id].ax,&mater[mat_id].ay,&mater[mat_id].az,
                &mater[mat_id].sig,&mater[mat_id].devPotFlag,&mater[mat_id].volPotFlag);
  }

  if(ferror(in)){
    PGFEM_printerr("ERROR: fscanf returned error in: %s(%s)\n",__func__,__FILE__);
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
    CHECK_SCANF(in,"%lf",&matgeom->cf[i]);
    matgeom->cm[i] = 1.0 - matgeom->cf[i];
    matgeom->cd[i] = 0.0;
  }
  for (i=0;i<np;i++){
    for (j=0;j<9;j++){
      CHECK_SCANF(in,"%lf",&matgeom->ee[i][j]);
    }
  }
  CHECK_SCANF(in,"%ld %lf %lf",&matgeom->SH,&matgeom->a1,&matgeom->a2);

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

    CHECK_SCANF(in,"%ld",&znod[i].nod);

    for (j=0;j<ndofn;j++)
      CHECK_SCANF(in,"%lf",&znod[i].load[j]);
  }

}

void read_elem_surface_load (FILE *in,
                             long nle_s,
                             long ndofn,
                             Element *elem,
                             ZATELEM *zele_s)
/*

 */
{
  long i,j,ii,jj{};

  for (i=0;i<nle_s;i++){
    CHECK_SCANF(in,"%ld",&zele_s[i].elem);
    ii = elem[zele_s[i].elem].toe;
    if (ii == 4)  jj = 3;
    if (ii == 8)  jj = 4;
    if (ii == 10) jj = 6;
    zele_s[i].sur = PGFEM_calloc (long, jj);
    for (j=0;j<jj;j++){
      CHECK_SCANF(in,"%ld",&zele_s[i].sur[j]);
    }
    for (j=0;j<ndofn;j++){
      CHECK_SCANF(in,"%lf",&zele_s[i].load[j]);
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
    CHECK_SCANF(in,"%lf",&sup->defl_d[i]);
  }
  err = ferror(in);
  PGFEM_fclose(in);
  return err;
}

int override_material_properties(const size_t nmat,
                                 const PGFem3D_opt *opt,
                                 Material *mater)
{
  int err = 0;
  int n_override = 0;
  int idx = -1;

  /* exit early if no override file is specified */
  if ( opt->override_material_props == NULL ) return err;

  /* else */
  FILE *in = PGFEM_fopen(opt->override_material_props,"r");
  err += scan_for_valid_line(in);
  if (err) goto exit_err;

  /* number of materials to override */
  CHECK_SCANF(in,"%d",&n_override);
  assert( n_override >= 0 );
  assert( (size_t)n_override <= nmat );

  for (int i=0; i<n_override; i++){
    /* scan for override data */
    err += scan_for_valid_line(in);
    if (err) goto exit_err;

    /* read override data */
    idx = -1; /* poisoned value */
    CHECK_SCANF(in,"%d",&idx);
    assert(idx >= 0);
    assert( (size_t)idx < nmat );
    err += read_material(in,idx,mater,opt->legacy);
  }

 exit_err:
  if( ferror(in) ) err ++;
  PGFEM_fclose(in);
  return err;
}
