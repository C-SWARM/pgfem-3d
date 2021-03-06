#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "element.h"
#include "allocation.h"
#include "enumerations.h"
#include "utils.h"

using namespace pgfem3d;

void Element::init() {
  pr = -1;
  toe = -1;
  nod = nullptr;
  mat = nullptr;
  hom = nullptr;
  bnd_type = nullptr;
  bnd_id = nullptr;
  n_dofs = 0;
  G_dof_ids = nullptr;
  L_dof_ids = nullptr;
  n_be = 0;
  be_ids = nullptr;
  n_bub = 0;
  n_bub_dofs = 0;
  bub_dofs = nullptr;
  d_bub_dofs = nullptr;
  L = nullptr;
  LO = nullptr;
}

Element* build_elem (FILE *in,
                     const long ne,
                     const int analysis)
/*
  ne - Number of elements
  type - Type of finite element
*/
{
  Element *pom = PGFEM_calloc (Element, ne);

  for (long i = 0; i < ne; ++i) {
    /* initialize all variables for current element. */
    pom[i].init();

    /* Read elem types */
    CHECK_SCANF (in,"%ld",&pom[i].toe);

    // Verify toe makes sense.
    switch (analysis) {
     case STABILIZED:
     case MINI:
     case MINI_3F:
      if (pom[i].toe != 4) {
        PGFEM_printf("Incorrect element type in build_elem\n");
        PGFEM_Abort();
      }
      break;
     case DISP:
      //      if (pom[i].toe != 4){
      //  PGFEM_printf("This element type has not been tested"
      //         " with this analysis type (yet).\n");
      //  PGFEM_Abort();
      //      }
      break;
     default:
      /* do nothing */
      break;
    }

    pom[i].nod = new long[pom[i].toe]{};
    pom[i].mat = new long[3]{};
    // elem[i].mat[0] -> Matrix
    // elem[i].mat[1] -> Fiber
    // elem[i].mat[2] -> Homogeneous
    pom[i].hom = new long[2]{};

    /* bounding features */
    switch (pom[i].toe) {
     case 4: /* tets */
     case 10:
      pom[i].bnd_type = new int[4]{};
      pom[i].bnd_id = new int[4]{};
      break;
     case 8:
      pom[i].bnd_type = new int[6]{};
      pom[i].bnd_id = new int[6]{};
      break;
     default:
      PGFEM_fprintf(stderr,"ERROR: Unsupported element type! %s:%s:%d\n",
                    __func__,__FILE__,__LINE__);
      PGFEM_Abort();
      break;
    }

    switch (analysis) {
     case MINI:
     case MINI_3F:
      pom[i].n_bub_dofs = 3;
      pom[i].n_bub = 1;
      pom[i].bub_dofs = new double[pom[i].n_bub * pom[i].n_bub_dofs]{};
      pom[i].d_bub_dofs = new double[pom[i].n_bub * pom[i].n_bub_dofs]{};
      break;

     default:
      // pom[i].bub_dofs = new double[1]{};
      // pom[i].d_bub_dofs = new double[1]{};
      break;
    }

    /* GFEM */
    /* if (gem == 1) {
       pom[i].gtyp = pom[i].gfem = pom[i].gele = -1;
       pom[i].gx1 = new double[2]{};
       pom[i].gx2 = new double[2]{};
       pom[i].gx3 = new double[2]{});
       } */
  }/* for each element */

  return pom;
} /* build_elem() */

void Element::fini() {
  delete [] nod;
  delete [] mat;
  delete [] hom;
  delete [] bnd_type;
  delete [] bnd_id;
  delete [] G_dof_ids;
  delete [] L_dof_ids;
  delete [] be_ids;
  delete [] bub_dofs;
  delete [] d_bub_dofs;
  delete [] LO;
}

void destroy_elem(Element *elem, const long ne)
{
  for (long i = 0; i < ne; ++i) {
    elem[i].fini();
  }
  PGFEM_free(elem);
}/* destroy_elem() */

/*=== Declare static helper functions ===*/
static void read_assign_elem_material(FILE *in, Element& elem) {
  long *mat = elem.mat;
  long *hom = elem.hom;

  /* matrix || fiber || volume fraction (cf) || fibre orientation (psi) */
  CHECK_SCANF (in,"%ld %ld %ld %ld",
               mat,mat+1,hom,hom+1);
}

static int is_elem_supported(const Element& elem, const SUPP_1& sup) {
  for (long j = 0; j < sup.ndn; ++j) {
    for (int k = 0; k < elem.toe; ++k) {
      if (sup.lnpd[j] == elem.nod[k]) {
        return 1;
      }
    }
  }
  return 0;
}

static void read_tet_conn(FILE *in, Element& elem) {
  /* get pointers to element data */
  long     *nod = elem.nod;
  long      *pr = &elem.pr;
  int *bnd_type = elem.bnd_type;
  int   *bnd_id = elem.bnd_id;
  long    ftype = 0;
  long      fid = 0;

  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld "
               "%d %d %d %d %d %d %d %d",
               nod,nod+1,nod+2,nod+3,&ftype,&fid,pr,
               bnd_id,bnd_id+1,bnd_id+2,bnd_id+3,
               bnd_type,bnd_type+1,bnd_type+2,bnd_type+3);
}

static void read_qtet_conn(FILE *in, Element& elem) {
  /* get pointers to element data */
  long     *nod = elem.nod;
  long      *pr = &elem.pr;
  int *bnd_type = elem.bnd_type;
  int   *bnd_id = elem.bnd_id;
  long    ftype = 0;
  long      fid = 0;

  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld "
               "%ld %ld %ld  %d %d %d %d "
               "%d %d %d %d",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,nod+8,nod+9,
               &ftype,&fid,pr,bnd_id,bnd_id+1,bnd_id+2,bnd_id+3,
               bnd_type,bnd_type+1,bnd_type+2,bnd_type+3);
}

static void read_hex_conn(FILE *in, Element& elem) {
  /* get pointers to element data */
  long     *nod = elem.nod;
  long      *pr = &elem.pr;
  int *bnd_type = elem.bnd_type;
  int   *bnd_id = elem.bnd_id;
  long    ftype = 0;
  long      fid = 0;

  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld "
               "%d %d %d %d %d %d "
               "%d %d %d %d %d %d",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,&ftype,&fid,pr,
               bnd_id,bnd_id+1,bnd_id+2,bnd_id+3,bnd_id+4,bnd_id+5,
               bnd_type,bnd_type+1,bnd_type+2,bnd_type+3,bnd_type+4,
               bnd_type+5);
}

/*=== LEGACY FILE FORMAT ===*/
static void read_tet_conn_legacy(FILE *in, Element& elem) {
  /* get pointers to element data */
  long *nod = elem.nod;
  long  *pr = &elem.pr;

  CHECK_SCANF (in,"%ld %ld %ld %ld %ld",
               nod,nod+1,nod+2,nod+3,pr);
}

static void read_qtet_conn_legacy(FILE *in, Element& elem) {
  /* get pointers to element data */
  long *nod = elem.nod;
  long  *pr = &elem.pr;

  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,nod+8,nod+9,pr);
}

static void read_hex_conn_legacy(FILE *in, Element& elem) {
  /* get pointers to element data */
  long *nod = elem.nod;
  long  *pr = &elem.pr;

  CHECK_SCANF (in,"%ld %ld %ld %ld %ld %ld %ld %ld %ld",
               nod,nod+1,nod+2,nod+3,nod+4,nod+5,nod+6,nod+7,pr);
}

/*
  in   - Input file
  ne   - Number of elements
  nl   - Number of layers
  toe  - Array of types of element
  elem - Structure type of Element

  %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
*/
void read_elem (FILE *in,
                long ne,
                Element *elem,
                SUPP sup,
                const int legacy)
{
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);

  for (long i = 0; i < ne ; ++i) {
    int err = 0;

    /* get elem id */
    long id = 0;
    CHECK_SCANF (in,"%ld",&id);
    assert(0 <= id and id < ne);

    if (legacy) {
      switch (elem[i].toe) {
       case 4:  read_tet_conn_legacy (in, elem[id]); break;
       case 10: read_qtet_conn_legacy(in, elem[id]); break;
       case 8:  read_hex_conn_legacy (in, elem[id]); break;
       default:
        err = 1;
        break;
      }
    } else {
      switch (elem[i].toe) {
       case 4:  read_tet_conn (in, elem[id]); break;
       case 10: read_qtet_conn(in, elem[id]); break;
       case 8:  read_hex_conn (in, elem[id]); break;
       default:
        err = 1;
        break;
      }
    } /* not legacy format */

    if (ferror(in)) {
      PGFEM_printerr("[%d]ERROR:CHECK_SCANF returned error"
                     " reading element %ld!\n",err_rank,i);
      PGFEM_Abort();
    }

    if (feof(in)) {
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                     err_rank);
      PGFEM_Abort();
    }

    if (err) {
      PGFEM_printerr("[%d]ERROR: element %ld is of unrecognized type (%ld)!\n",
                     err_rank, i, elem[i].toe);
      PGFEM_Abort();
    }
  }

  for (long i = 0; i < ne; ++i) {
    long id;
    CHECK_SCANF (in, "%ld", &id);
    assert(0 <= id and id < ne);
    read_assign_elem_material(in, elem[id]);
    if (ferror(in) || feof(in)) {
      PGFEM_printerr("[%d]ERROR: CHECK_SCANF error assigning element "
                     "material!\n", err_rank);
      PGFEM_Abort();
    }

    /* check if element is supported */
    sup->nde += is_elem_supported(elem[i], *sup);
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


