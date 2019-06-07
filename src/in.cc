#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "allocation.h"
#include "in.h"
#include "utils.h"
#include <cstring>
#include <cassert>

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

using namespace pgfem3d;

/// read Dirichlet boundary conditions on nodes
///
/// This function will generate and return SUPP object.
/// if in==NULL, assumed no boundary conditions. Code will proceed to next
/// step in reading boundary condition values.
///
/// \param[in]  in    Input file
/// \param[in]  ndofn Number of degrees of freedom in one node
/// \param[in]  node  Structure type of Node
/// \param[in]  mp_id multiphysics id
/// \return     sup   created boundary condition structure
SUPP read_Dirichlet_BCs(FILE *in,
                        long nn,
                        long ndofn,
                        Node *node,
                        const int mp_id)
{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if (PFEM_DEBUG) {
    PGFEM_printf("[%d] reading supports.\n", err_rank);
  }

  // Return 1 if the node with id `n` has a prescribed deflection.
  auto has_prescribed_deflection = [ndofn,node,mp_id](auto n) {
    auto& ids = node[n].id_map[mp_id].id;
    for (int k = 0; k < ndofn; ++k) {
      if (ids[k] == 1 || ids[k] < 0) {
        return 1;
      }
    }
    return 0;
  };

  SUPP sup = PGFEM_calloc (SUPP_1, 1);

  // read supported nodes
  if (in) {
    CHECK_SCANF(in, "%ld", &sup->nsn);
  }

  if (sup->nsn == 0)  {
    sup->supp = PGFEM_calloc (long, 1);
  }
  else {
    sup->supp = PGFEM_calloc (long, sup->nsn);
  }

  for (int i = 0, e = sup->nsn; i < e; ++i) {
    long n;
    CHECK_SCANF(in,"%ld", &n);
    sup->supp[i] = n;
    for (int k = 0; k < ndofn; ++k) {
      CHECK_SCANF(in, "%ld", &node[n].id_map[mp_id].id[k]);
    }
    sup->ndn += has_prescribed_deflection(n);

    if (ferror(in)) {
      PGFEM_printerr("[%d]ERROR:fscanf returned error reading support %ld!\n",
                     err_rank, i);
      PGFEM_Abort();
    }

    if (feof(in)) {
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                     err_rank);
      PGFEM_Abort();
    }
  }

  /***************************************************/
  /* create list of nodes with prescribed deflection */
  /***************************************************/
  if (sup->ndn == 0) {
    sup->lnpd = PGFEM_calloc (long, 1);
  }
  else {
    sup->lnpd = PGFEM_calloc (long, sup->ndn);
  }

  // Copy the list of nodes with prescribed deflection from `sup->supp` to
  // `sup->lnpd`.
  for (int i = 0, insert = 0, e = sup->nsn; i < e; ++i) {
    long n = sup->supp[i];
    if (has_prescribed_deflection(n)) {
      assert(insert < sup->ndn);
      sup->lnpd[insert++] = n;
    }
  }

  // allocate the prescribed macro deformation gradient
  sup->F0 = PGFEM_calloc(double, 9);
  sup->N0 = PGFEM_calloc(double, 3);

  return sup;
}

/// read Dirichlet boundary condition values
///
/// Before read BC values, SUPP object should be created.
/// Prior to call this function, read_Dirichlet_BCs function should be called first.
///
/// \param[in]  in    Input file
/// \param[in]  ndofn Number of degrees of freedom in one node
/// \param[in]  node  Structure type of Node
/// \param[in]  mp_id multiphysics id
/// \return non-zero on internal ERROR
int read_Dirichlet_BCs_values(FILE *in,
                              long nn,
                              long ndofn,
                              Node *node,
                              SUPP sup,
                              const int mp_id)
{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if (PFEM_DEBUG) {
    PGFEM_printf("[%d] reading BCs_values.\n", err_rank);
  }

  // read nodes with prescribed deflection
  CHECK_SCANF(in, "%ld", &sup->npd);
  assert(0 <= sup->npd);

  if (sup->npd == 0) {
    sup->defl   = PGFEM_calloc (double, 1);
    sup->defl_d = PGFEM_calloc (double, 1);
  }
  else {
    sup->defl   = PGFEM_calloc (double, sup->npd);
    sup->defl_d = PGFEM_calloc (double, sup->npd);
  }

  for (long i = 0; i < sup->npd; ++i) {
    CHECK_SCANF(in, "%lf", &sup->defl_d[i]);
  }

  if (ferror(in)) {
    PGFEM_printerr("[%d]ERROR:fscanf returned error"
                   " reading prescribed deflections!\n",err_rank);
    PGFEM_Abort();
  }
  if (feof(in)) {
    PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                   err_rank);
    PGFEM_Abort();
  }

  return 0;
}

SUPP read_supports(FILE *in,
                   long nn,
                   long ndofn,
                   Node *node,
                   const int mp_id)
/*
  in    - Input file
  nl    - Number of layers
  ndofn - Number of degrees of freedom in one node
  node  - Structure type of Node
  sup   - Structure type of SUPP

  %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
*/
{
  SUPP sup = read_Dirichlet_BCs(in, nn, ndofn, node, mp_id);
  if (read_Dirichlet_BCs_values(in, nn, ndofn, node, sup, mp_id)) {
    PGFEM_printf("Error reading Dirichlet boundary conditions");
    PGFEM_Abort();
  }
  return sup;
}



int read_material (FILE *in,
                   Material& mater,
                   const int legacy)
{
  if (legacy) {
    mater.devPotFlag = 0;
    mater.volPotFlag = 0;
    CHECK_SCANF(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &mater.Ex,  &mater.Ey,  &mater.Ez,
                &mater.Gyz, &mater.Gxz, &mater.Gxy,
                &mater.nyz, &mater.nxz, &mater.nxy,
                &mater.ax,  &mater.ay,  &mater.az,
                &mater.sig);
  } else {
    CHECK_SCANF(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
                &mater.Ex,  &mater.Ey,         &mater.Ez,
                &mater.Gyz, &mater.Gxz,        &mater.Gxy,
                &mater.nyz, &mater.nxz,        &mater.nxy,
                &mater.ax,  &mater.ay,         &mater.az,
                &mater.sig, &mater.devPotFlag, &mater.volPotFlag);
  }

  if (ferror(in)) {
    PGFEM_printerr("ERROR: fscanf returned error in: %s(%s)\n",__func__,__FILE__);
    return 1;
  }

  return 0;
}

void read_matgeom (FILE *in,
                   long nc,
                   long np,
                   MATGEOM matgeom)
{
  if (PFEM_DEBUG) {
    PGFEM_printf("[%d] reading material geom.\n", 1);
  }

  for (long i = 0; i < nc; ++i) {
    CHECK_SCANF(in,"%lf", &matgeom->cf[i]);
    matgeom->cm[i] = 1.0 - matgeom->cf[i];
    matgeom->cd[i] = 0.0;
  }
  for (long i = 0; i < np; ++i) {
    for (int j = 0; j < 9; ++j) {
      CHECK_SCANF(in,"%lf", &matgeom->ee[i][j]);
    }
  }
  CHECK_SCANF(in,"%ld %lf %lf", &matgeom->SH, &matgeom->a1, &matgeom->a2);
}

void read_nodal_load (FILE *in,
                      long nln,
                      long ndofn,
                      ZATNODE *znod)
{
  if (PFEM_DEBUG) {
    PGFEM_printf("[%d] reading nodal loads.\n",1);
  }

  for (long i = 0; i < nln; ++i) {
    CHECK_SCANF(in, "%ld", &znod[i].nod);
    for (long j = 0; j < ndofn; ++j) {
      CHECK_SCANF(in, "%lf", &znod[i].load[j]);
    }
  }
}

void read_elem_surface_load (FILE *in,
                             long nle_s,
                             long ndofn,
                             Element *elem,
                             ZATELEM *zele_s)
{
  for (long i = 0; i < nle_s; ++i) {
    CHECK_SCANF(in, "%ld", &zele_s[i].elem);
    long jj;
    switch (elem[zele_s[i].elem].toe) {
     case 4:  jj = 3; break;
     case 8:  jj = 4; break;
     case 10: jj = 6; break;
     default:
      PGFEM_printerr("Invalid type of element\n");
      PGFEM_Abort();
    }

    zele_s[i].sur = PGFEM_calloc (long, jj);
    for (long j = 0; j < jj; ++j) {
      CHECK_SCANF(in, "%ld", &zele_s[i].sur[j]);
    }
    for (long j = 0; j < ndofn; ++j) {
      CHECK_SCANF(in, "%lf", &zele_s[i].load[j]);
    }
  }

}

int override_prescribed_displacements(SUPP_1& sup, const char* fn) {
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  if (err_rank == 0) {
    PGFEM_printf("Overriding the prescribed displacements with:\n%s\n", fn);
  }

  if (auto in = scoped_fopen(fn, "r")) {
    for (long i = 0, e = sup.npd; i < e; ++i) {
      CHECK_SCANF(in, "%lf", &sup.defl_d[i]);
    }
    return 0;
  }
  return 1;
}

int override_material_properties(const long nmat,
                                 const PGFem3D_opt *opt,
                                 Material *mater)
{
  /* exit early if no override file is specified */
  if (opt->override_material_props == nullptr) {
    return 0;
  }

  FILE *in = PGFEM_fopen(opt->override_material_props,"r");
  if (in == nullptr) {
    std::cerr << "Could not open " << opt->override_material_props << "\n";
    return 1;
  }

  /* number of materials to override */
  int n_override = 0;
  if (scan_for_valid_line(in)) {
    goto exit_err;
  }
  CHECK_SCANF(in, "%d", &n_override);
  assert(0 <= n_override and n_override < nmat);

  for (int i = 0, e = n_override; i < e; ++i) {
    if (scan_for_valid_line(in)) {
      goto exit_err;
    }

    int idx = -1;
    CHECK_SCANF(in, "%d", &idx);
    assert(0 <= idx and idx < nmat);
    if (read_material(in, mater[idx], opt->legacy)) {
      goto exit_err;
    }
  }

 exit_err:
  int err = (ferror(in)) ? 2 : 1;
  PGFEM_fclose(in);
  return err;
}
