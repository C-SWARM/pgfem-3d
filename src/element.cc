#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "element.h"
#include "allocation.h"
#include "enumerations.h"
#include "utils.h"

using namespace pgfem3d;

Element* build_elem (FILE *in,
                     const long ne,
                     const int analysis)
/*
  ne - Number of elements
  type - Type of finite element
*/
{
  Element *pom = PGFEM_calloc (Element, ne);

  for (long ii=0;ii<ne;ii++){
    /* initialize all variables for current element. */
    pom[ii].pr = -1;
    pom[ii].toe = -1;
    pom[ii].nod = NULL;
    pom[ii].mat = NULL;
    pom[ii].hom = NULL;
    pom[ii].bnd_type = NULL;
    pom[ii].bnd_id = NULL;
    pom[ii].n_dofs = 0;
    pom[ii].G_dof_ids = NULL;
    pom[ii].L_dof_ids = NULL;
    pom[ii].n_be = 0;
    pom[ii].be_ids = NULL;
    pom[ii].n_bub = 0;
    pom[ii].n_bub_dofs = 0;
    pom[ii].bub_dofs = NULL;
    pom[ii].d_bub_dofs = NULL;
    pom[ii].L = NULL;
    pom[ii].LO = NULL;

    /* Read elem types */
    CHECK_SCANF (in,"%ld",&pom[ii].toe);

    switch(analysis){
     case STABILIZED:
     case MINI:
     case MINI_3F:
      if (pom[ii].toe != 4){
        PGFEM_printf("Incorrect element type in  build_elem\n");
        PGFEM_Abort();
      }
      break;
     case DISP:
      //      if (pom[ii].toe != 4){
      //  PGFEM_printf("This element type has not been tested"
      //         " with this analysis type (yet).\n");
      //  PGFEM_Abort();
      //      }
      break;
     default:
      /* do nothing */
      break;
    }

    pom[ii].nod = PGFEM_calloc (long, pom[ii].toe);
    pom[ii].mat = PGFEM_calloc (long, 3);
    // elem[ii].mat[0] -> Matrix
    // elem[ii].mat[1] -> Fiber
    // elem[ii].mat[2] -> Homogeneous
    pom[ii].hom = PGFEM_calloc (long, 2);

    /* bounding features */
    switch(pom[ii].toe){
     case 4: /* tets */
     case 10:
      pom[ii].bnd_type = PGFEM_calloc(int, 4);
      pom[ii].bnd_id = PGFEM_calloc(int, 4);
      break;
     case 8:
      pom[ii].bnd_type = PGFEM_calloc(int, 6);
      pom[ii].bnd_id = PGFEM_calloc(int, 6);
      break;
     default:
      PGFEM_fprintf(stderr,"ERROR: Unsupported element type! %s:%s:%d\n",
                    __func__,__FILE__,__LINE__);
      PGFEM_Abort();
      break;
    }

    switch(analysis){
     case MINI:
     case MINI_3F:
      pom[ii].n_bub_dofs = 3;
      pom[ii].n_bub = 1;
      pom[ii].bub_dofs = PGFEM_calloc(double, pom[ii].n_bub*pom[ii].n_bub_dofs);
      pom[ii].d_bub_dofs = PGFEM_calloc(double, pom[ii].n_bub*pom[ii].n_bub_dofs);
      break;

     default:
      pom[ii].bub_dofs = PGFEM_calloc (double, 1);
      pom[ii].d_bub_dofs = PGFEM_calloc (double, 1);
      break;
    }

    /* GFEM */
    /* if (gem == 1) {
       pom[ii].gtyp = pom[ii].gfem = pom[ii].gele = -1;
       pom[ii].gx1 = PGFEM_calloc (double, 2);
       pom[ii].gx2 = PGFEM_calloc (double, 2);
       pom[ii].gx3 = PGFEM_calloc (double, 2);
       } */
  }/* for each element */

  return (pom);
} /* build_elem() */

void destroy_elem(Element *elem,
                  const long ne)
{
  for(long i=0; i<ne; i++){
    free(elem[i].nod);
    free(elem[i].mat);
    free(elem[i].hom);
    free(elem[i].bnd_type);
    free(elem[i].bnd_id);
    free(elem[i].G_dof_ids);
    free(elem[i].L_dof_ids);
    free(elem[i].be_ids);
    free(elem[i].bub_dofs);
    free(elem[i].d_bub_dofs);
    free(elem[i].LO);
  }
  free(elem);
}/* destroy_elem() */

/*=== Declare static read functions ===*/
static int read_tet_conn_legacy(FILE *in,Element *elem);
static int read_qtet_conn_legacy(FILE *in,Element *elem);
static int read_hex_conn_legacy(FILE *in,Element *elem);
static int read_tet_conn(FILE *in,Element *elem);
static int read_qtet_conn(FILE *in,Element *elem);
static int read_hex_conn(FILE *in,Element *elem);
static int read_assign_elem_material(FILE *in,Element *elem);
static int is_elem_supported(const Element *p_elem,
                             const long nsup_node,
                             const long *sup_node_id);

void read_elem (FILE *in,
                long ne,
                Element *elem,
                SUPP sup,
                const int legacy)
/*
  in   - Input file
  ne   - Number of elements
  nl   - Number of layers
  toe  - Array of types of element
  elem - Structure type of Element

  %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
*/

{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);

  int err = 0;

  for (long i=0;i<ne;i++){
    if(legacy){
      switch(elem[i].toe){
       case 4: read_tet_conn_legacy(in,elem); break;
       case 10: read_qtet_conn_legacy(in,elem); break;
       case 8: read_hex_conn_legacy(in,elem); break;
       default: err = 1; break;
      }
    } else {
      switch(elem[i].toe){
       case 4: read_tet_conn(in,elem); break;
       case 10: read_qtet_conn(in,elem); break;
       case 8: read_hex_conn(in,elem); break;
       default: err = 1; break;
      }
    } /* not legacy format */

    if(ferror(in)){
      PGFEM_printerr("[%d]ERROR:CHECK_SCANF returned error"
                     " reading element %ld!\n",err_rank,i);
      PGFEM_Abort();
    } else if(feof(in)){
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                     err_rank);
      PGFEM_Abort();
    } else if(err){
      PGFEM_printerr("[%d]ERROR: element %ld is of"
                     " unrecognized type (%ld)!\n",err_rank,i,elem[i].toe);
      PGFEM_Abort();
    }
  }

  for (long i=0;i<ne;i++){
    err = read_assign_elem_material(in,elem);
    if(err){
      PGFEM_printerr("[%d]ERROR: CHECK_SCANF error assigning "
                     "element material!\n",err_rank);
      PGFEM_Abort();
    }

    /* check if element is supproted */
    sup->nde += is_elem_supported(&elem[i],sup->ndn,sup->lnpd);
  }/* end i<ne */
}/* read_elem() */

void write_element_fname(const char *filename,
                         const int nelem,
                         const Element *elems)
{
  FILE *ofile = fopen(filename,"w");
  if(ofile == NULL){
    PGFEM_printerr("ERROR: cannot open file %s in %s\n",filename,__func__);
    PGFEM_Abort();
  }

  write_element(ofile,nelem,elems);

  fclose(ofile);
} /* write_elem_fname() */

void write_element(FILE *ofile,
                   const int nelem,
                   const Element *elems)
{
  /* print header describing format */
  PGFEM_fprintf(ofile,"Each element is printed in a block with each line in\n"
                "the block formatted as follows:\n"
                "   id nne   nod\n"
                "n_dof       Lid::Gid ...\n"
                " n_be       be_ids...\n");
  PGFEM_fprintf(ofile,"====================================================\n");
  for(int i=0; i<nelem; i++){
    const Element *p_elem = &elems[i];

    /* line 1 (connectivity */
    PGFEM_fprintf(ofile,"%5d %3ld ",i,p_elem->toe);
    for(int j=0; j<p_elem->toe; j++){
      PGFEM_fprintf(ofile,"%5ld ",p_elem->nod[j]);
    }
    PGFEM_fprintf(ofile,"\n");

    /* line 2 (elem dofs) */
    PGFEM_fprintf(ofile,"%5d     ",p_elem->n_dofs);
    for(int j=0; j<p_elem->n_dofs; j++){
      PGFEM_fprintf(ofile,"%5ld::%-5ld ",p_elem->L_dof_ids[j],p_elem->G_dof_ids[j]);
    }
    PGFEM_fprintf(ofile,"\n");

    /* line 3 (bounding elements) */
    PGFEM_fprintf(ofile,"%5d     ",p_elem->n_be);
    for(int j=0; j<p_elem->n_be; j++){
      PGFEM_fprintf(ofile,"%5ld ",p_elem->be_ids[j]);
    }
    PGFEM_fprintf(ofile,"\n\n");
  }
}

static int read_assign_elem_material(FILE *in,
                                     Element *elem)
{
  int err = 0;
  long id;
  CHECK_SCANF (in,"%ld",&id);
  long *mat = elem[id].mat;
  long *hom = elem[id].hom;

  /* matrix || fiber || volume fraction (cf) || fibre orientation (psi) */
  CHECK_SCANF (in,"%ld %ld %ld %ld",mat,mat+1,hom,hom+1);

  if(ferror(in) || feof(in)) err = 1;
  return err;
}/* read_assign_elem_material */

static int is_elem_supported(const Element *p_elem,
                             const long nsup_node,
                             const long *sup_node_id)
{
  const long *nod = p_elem->nod;
  const int nne = p_elem->toe;

  for (long j=0; j<nsup_node; j++){
    for (int k=0; k<nne; k++){
      if (sup_node_id[j] == nod[k]){
        return 1;
      }
    }
  }

  return 0;
}

static int read_tet_conn(FILE *in,
                         Element *elem)
{
  int err = 0;

  /* get elem id */
  long id = 0;
  CHECK_SCANF (in,"%ld",&id);

  /* get pointers to element data */
  long *nod = elem[id].nod;
  long *pr = &elem[id].pr;
  int *bnd_type = elem[id].bnd_type;
  int *bnd_id = elem[id].bnd_id;
  long ftype = 0;
  long fid = 0;

  /* read data */
  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld "
               "%d %d %d %d %d %d %d %d",
               nod,nod+1,nod+2,nod+3,&ftype,&fid,pr,
               bnd_id,bnd_id+1,bnd_id+2,bnd_id+3,
               bnd_type,bnd_type+1,bnd_type+2,bnd_type+3);

  /* check for file error */
  if(ferror(in) != 0) err = 1;

  return err;
}/* read_tet_conn() */

static int read_qtet_conn(FILE *in,
                          Element *elem)
{
  int err = 0;

  /* get elem id */
  long id = 0;
  CHECK_SCANF (in,"%ld",&id);

  /* get pointers to element data */
  long *nod = elem[id].nod;
  long *pr = &elem[id].pr;
  int *bnd_type = elem[id].bnd_type;
  int *bnd_id = elem[id].bnd_id;
  long ftype = 0;
  long fid = 0;

  /* read data */
  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld "
               "%ld %ld %ld  %d %d %d %d "
               "%d %d %d %d",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,nod+8,nod+9,
               &ftype,&fid,pr,bnd_id,bnd_id+1,bnd_id+2,bnd_id+3,
               bnd_type,bnd_type+1,bnd_type+2,bnd_type+3);

  /* check for file error */
  if(ferror(in) != 0) err = 1;

  return err;
}/* read_qtet_conn() */

static int read_hex_conn(FILE *in,
                         Element *elem)
{
  int err = 0;

  /* get elem id */
  long id = 0;
  CHECK_SCANF (in,"%ld",&id);

  /* get pointers to element data */
  long *nod = elem[id].nod;
  long *pr = &elem[id].pr;
  int *bnd_type = elem[id].bnd_type;
  int *bnd_id = elem[id].bnd_id;
  long ftype = 0;
  long fid = 0;

  /* read data */
  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld "
               "%d %d %d %d %d %d "
               "%d %d %d %d %d %d",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,&ftype,&fid,pr,
               bnd_id,bnd_id+1,bnd_id+2,bnd_id+3,bnd_id+4,bnd_id+5,
               bnd_type,bnd_type+1,bnd_type+2,bnd_type+3,bnd_type+4,
               bnd_type+5);

  /* check for file error */
  if(ferror(in) != 0) err = 1;

  return err;
}/* read_hex_conn() */


/*=== LEGACY FILE FORMAT ===*/
static int read_tet_conn_legacy(FILE *in,
                                Element *elem)
{
  int err = 0;

  /* get elem id */
  long id = 0;
  CHECK_SCANF (in,"%ld",&id);

  /* get pointers to element data */
  long *nod = elem[id].nod;
  long *pr = &elem[id].pr;

  /* read data */
  CHECK_SCANF (in,"%ld %ld %ld %ld %ld",nod,nod+1,nod+2,nod+3,pr);

  /* check for file error */
  if(ferror(in) != 0) err = 1;

  return err;
}/* read_tet_conn_legacy() */

static int read_qtet_conn_legacy(FILE *in,
                                 Element *elem)
{
  int err = 0;

  /* get elem id */
  long id = 0;
  CHECK_SCANF (in,"%ld",&id);

  /* get pointers to element data */
  long *nod = elem[id].nod;
  long *pr = &elem[id].pr;

  /* read data */
  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,nod+8,nod+9,pr);

  /* check for file error */
  if(ferror(in) != 0) err = 1;

  return err;
}/* read_qtet_conn_legacy() */

static int read_hex_conn_legacy(FILE *in,
                                Element *elem)
{
  int err = 0;

  /* get elem id */
  long id = 0;
  CHECK_SCANF (in,"%ld",&id);

  /* get pointers to element data */
  long *nod = elem[id].nod;
  long *pr = &elem[id].pr;

  /* read data */
  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,pr);

  /* check for file error */
  if(ferror(in) != 0) err = 1;

  return err;
}/* read_hex_conn_legacy() */

