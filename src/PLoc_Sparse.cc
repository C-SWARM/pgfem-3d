/* HEADER */
#include "PLoc_Sparse.h"
#include "allocation.h"
#include "enumerations.h"
#include "utils.h"
#include <string.h>

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

using pgfem3d::SparseComm;
using pgfem3d::solvers::SparseSystem;

void PLoc_Sparse (double **Lk,
                  double *lk,
                  Ai_t *Ai,
                  int *Ap,
                  long *cnL,
                  long *cnG,
                  long ndofe,      // should be int
                  long *Ddof,
                  long GDof,
                  int myrank,
                  int nproc,
                  SparseComm *comm,
                  int interior,
                  SparseSystem *system,
                  const int analysis)
{
  int err = 0;

  int nrows = 0;
  int ncols = 0;
  int nsend = 0;

  SparseSystem::Index *gDofID; 
  int *send, *rows, *cols;

  /* The element stiffness matrix is by definition ndofe x ndofe,
     therefore we start with that many rows and columns and subtract
     rows if there is a BC (gDofID < 0) or if the dof is in another
     domain.  Columns are only removed if there is a BC on the dof. */

  gDofID = aloc1l(ndofe);
  send = aloc1i(ndofe);
  rows = aloc1i(ndofe);
  cols = aloc1i(ndofe);


  for(int i = 0; i < (int) ndofe; ++i){
    gDofID[i] = cnG[i] - 1;                 // varies bewtween: negative < cnG < long 
    if(gDofID[i] < 0){
      continue;
    }
    if (system->isLocal(gDofID[i])) {
      rows[nrows++] = i;
      cols[ncols++] = i;
    } else {
      send[nsend++] = i;
      cols[ncols++] = i;
    }
  }

  if(interior && nsend){
    PGFEM_printf("[%d]WARNING: interior element needs to send %d rows?!?",
                 myrank,nsend);
  }

  /* Allocate the rows, columns, and values */
  SparseSystem::Index *row_idx, *n_cols, *col_idx;
  double *values;
  if(nrows > 0){
    row_idx = PGFEM_calloc(SparseSystem::Index, nrows);
    n_cols = PGFEM_calloc(SparseSystem::Index, nrows);
  } else {
    row_idx = PGFEM_calloc(SparseSystem::Index, 1);
    n_cols = PGFEM_calloc(SparseSystem::Index, 1);
  }

  if((ncols > 0) && (nrows >0)){
    col_idx = PGFEM_calloc(SparseSystem::Index, nrows*ncols);
    values = aloc1(nrows*ncols);
  } else {
    col_idx = PGFEM_calloc(SparseSystem::Index, 1);
    values = aloc1(1);
  }

  for(int i = 0; i < nrows; ++i){
    row_idx[i] = gDofID[rows[i]];
    n_cols[i] = ncols;
    for(int j = 0; j < ncols; ++j){
      values[i*ncols + j] = lk[rows[i]*ndofe + cols[j]];
      col_idx[i*ncols + j] = gDofID[cols[j]];
    }
  }

  if(nrows > 0){
    if(PFEM_DEBUG){
      char ofile[50];
      switch(analysis){
       case STABILIZED:
        sprintf(ofile,"stab_assem_loc_%d.log",myrank);
        break;
       case MINI:
        sprintf(ofile,"MINI_assem_loc_%d.log",myrank);
        break;
       case MINI_3F:
        sprintf(ofile,"MINI_3f_assem_loc_%d.log",myrank);
        break;
       default:
        sprintf(ofile,"el_assem_loc_%d.log",myrank);
        break;
      }

      FILE *out;
      out = fopen(ofile,"a");
      PGFEM_fprintf(out,"********************************************\n");
      print_array_i(out,n_cols,nrows,1,nrows);
      print_array_i(out,row_idx,nrows,1,nrows);
      print_array_i(out,col_idx,nrows*ncols,1,nrows*ncols);
      print_array_d(out,values,nrows*ncols,1,nrows*ncols);
      fclose(out);
    }
    /* Add to the stiffness matrix */
    system->add(nrows, n_cols, row_idx, col_idx, values);
  }

  if(err){
    PGFEM_printerr("[%d] Recieved error from HYPRE_IJMatrixAddToValues\n",myrank);
    print_array_i(PGFEM_stderr,n_cols,nrows,1,nrows);
    print_array_i(PGFEM_stderr,row_idx,nrows,1,nrows);
    print_array_i(PGFEM_stderr,col_idx,nrows*ncols,nrows,ncols);
  }

  /* free unnecessary arrays */
  PGFEM_free(row_idx);
  PGFEM_free(col_idx);
  PGFEM_free(values);
  PGFEM_free(n_cols);
  PGFEM_free(rows);

  /* Add to the send portion */
  for(int i=0; i<nsend; i++){
    long gID = gDofID[send[i]];
    long proc = 0;
    long srow_id = 0;
    long s_nnz = 0;
    for(long j=0; j<comm->Ns; j++){
      proc = comm->Nss[j];
      for(srow_id=0; srow_id<comm->S[proc]; srow_id++){
        if(gID == comm->LG[comm->SLID[proc][srow_id]]){
          goto found_matching_gID;
        }
      }
    }
   found_matching_gID: /* label :: breaks two for-loops when matching
                          gID found */
    for(long j=0; j<srow_id; j++){
      s_nnz += comm->SAp[proc][j];
    }
    for(int j = 0; j < ncols; ++j){
      gID = gDofID[cols[j]];
      for(long k=0; k<comm->SAp[proc][srow_id]; k++){
        if(gID == comm->SGRId[proc][s_nnz+k]){
          Lk[proc][s_nnz+k] += lk[send[i]*ndofe +cols[j]];
          break;
        }
      }
    }
  }

  /* free remaining memory */
  PGFEM_free(cols);
  PGFEM_free(send);
  PGFEM_free(gDofID);
}
