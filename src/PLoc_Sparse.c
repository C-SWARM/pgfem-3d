/* HEADER */
#include "PLoc_Sparse.h"
#include <string.h>
#include "PGFEM_mpi.h"
#include "enumerations.h"
#include "allocation.h"
#include "utils.h"

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

void PLoc_Sparse (double **Lk,
		  double *lk,
		  int *Ai,
		  int *Ap,
		  long *cnL,
		  long *cnG,
		  long ndofe,
		  int *Ddof,
		  long GDof,
		  int myrank,
		  int nproc,
		  COMMUN comm,
		  int interior,
		  PGFEM_HYPRE_solve_info *PGFEM_hypre,
		  const int analysis)
{
  long LI1 = PGFEM_hypre->ilower;
  long LI2 = PGFEM_hypre->iupper;

  int err = 0;

  int nrows = 0;
  int ncols = 0;
  int nsend = 0;


  int *gDofID, *send, *rows, *cols;

  /* The element stiffness matrix is by definition ndofe x ndofe,
     therefore we start with that many rows and columns and subtract
     rows if there is a BC (gDofID < 0) or if the dof is in another
     domain.  Columns are only removed if there is a BC on the dof. */

  gDofID = aloc1i(ndofe);
  send = aloc1i(ndofe);
  rows = aloc1i(ndofe);
  cols = aloc1i(ndofe);

  for(int i=0; i< (int) ndofe;i++){
    gDofID[i] = (int) (cnG[i] - 1);
    if(gDofID[i] < 0){
      continue;
    }
    if((LI1 <= gDofID[i]) && (gDofID[i] <= LI2)){
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
  int *row_idx, *n_cols, *col_idx;
  double *values;
  if(nrows > 0){
    row_idx = aloc1i(nrows);
    n_cols = aloc1i(nrows);
  } else {
    row_idx = aloc1i(1);
    n_cols = aloc1i(1);
  }

  if((ncols > 0) && (nrows >0)){
    col_idx = aloc1i(nrows*ncols);
    values = aloc1(nrows*ncols);
  } else {
    col_idx = aloc1i(1);
    values = aloc1(1);
  }

  for(int i=0; i<nrows; i++){
    row_idx[i] = gDofID[rows[i]];
    n_cols[i] = ncols;
    for(int j=0; j<ncols; j++){
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
   err = HYPRE_IJMatrixAddToValues(PGFEM_hypre->hypre_k,nrows,n_cols,
				   row_idx,col_idx,values);
  }

  if(err){
    PGFEM_printerr("[%d] Recieved error from HYPRE_IJMatrixAddToValues\n",myrank);
    print_array_i(PGFEM_stderr,n_cols,nrows,1,nrows);
    print_array_i(PGFEM_stderr,row_idx,nrows,1,nrows);
    print_array_i(PGFEM_stderr,col_idx,nrows*ncols,nrows,ncols);
  }

  /* free unnecessary arrays */
  free(row_idx);
  free(col_idx);
  free(values);
  free(n_cols);
  free(rows);

  /* Add to the send portion */
  for(int i=0; i<nsend; i++){
    long gID = (long) gDofID[send[i]];
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
    for(int j=0; j<ncols; j++){
      gID = (long) gDofID[cols[j]];
      for(long k=0; k<comm->SAp[proc][srow_id]; k++){
	if(gID == comm->SGRId[proc][s_nnz+k]){
	  Lk[proc][s_nnz+k] += lk[send[i]*ndofe +cols[j]];
	  break;
	}
      }
    }
  }

  /* free remaining memory */
  free(cols);
  free(send);
  free(gDofID);
}

/** PLoc_Sparse_rec for general (non-square) matrices. */
void PLoc_Sparse_rec (double **Lk,
		      double *lk,
		      int *Ai,
		      int *Ap,
		      long *cnG_row,
		      long *cnG_col,
		      long nrow,
		      long ncol,
		      int *Ddof,
		      long GDof,
		      int myrank,
		      int nproc,
		      COMMUN comm,
		      int interior,
		      PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  /* This function is much the same as PLoc_Sparse with the additional
     ability of handling non-square matrices. */

  /* ilower & iupper are global variables from hypre_global.h */
  long LI1 = PGFEM_hypre->ilower;
  long LI2 = PGFEM_hypre->iupper;

  int nrows = 0;
  int ncols = 0;
  int nsend = 0;

  /* We subtract rows if the DOF is owned by another domain or is
     prescribed by a boundary condition. Columns are only removed if
     the DOF is prescribed by a BC. */

  int *G_row_id = aloc1i(nrow);
  int *G_col_id = aloc1i(ncol);
  int *send = aloc1i(nrow);
  int *rows = aloc1i(nrow);
  int *cols = aloc1i(ncol);

  /* loop through rows and increment number to keep/send */
  for(int i=0; i< nrow; i++){
    G_row_id[i] = (cnG_row[i] - 1);
    if(G_row_id[i] < 0) continue; /* BC */
    else if( (LI1 <= G_row_id[i]) && (G_row_id[i] <= LI2) ){
      rows[nrows++] = i;
    } else {
      send[nsend++] = i;
    }
  }

  /* loop through columns and keep any without BC */
  for(int i=0; i<ncol; i++){
    G_col_id[i] = (cnG_col[i] - 1);
    if(G_col_id[i] < 0) continue; /* BC */
    else cols[ncols++] = i;
  }

  /* now we know how many rows/columns to keep/send and what index
     they belong to. Now allocate space for the row and column values
     noting that the local matrix is dense. */
  int *row_idx, *n_cols, *col_idx;
  double *values;
  if(nrows > 0){
    row_idx = aloc1i(nrows);
    n_cols = aloc1i(nrows);
  } else {
    row_idx = aloc1i(1);
    n_cols = aloc1i(1);
  }

  if( (ncols > 0) && (nrows > 0)){
    col_idx = aloc1i(nrows*ncols);
    values = aloc1(nrows*ncols);
  } else {
    col_idx = aloc1i(1);
    values = aloc1(1);
  }
 
  /* fill local sparse matrix structure */
  for(int i=0; i<nrows; i++){
    row_idx[i] = G_row_id[rows[i]];
    n_cols[i] = ncols;
    for(int j=0; j<ncols; j++){
      int idx = i*ncols+j;
      col_idx[idx] = G_col_id[cols[j]];
      values[idx] = lk[rows[i]*ncols + cols[j]];
    }
  }

  if(nrows > 0){
    HYPRE_IJMatrixAddToValues(PGFEM_hypre->hypre_k,nrows,n_cols,
			      row_idx,col_idx,values);
  }

  /* free unnecessary arrays */
  free(row_idx);
  free(col_idx);
  free(values);
  free(n_cols);
  free(rows);

  /* Add to the send portion */
  for(int i=0; i<nsend; i++){
    long gID = G_row_id[send[i]];
    long proc = 0;
    long srow_id = 0;
    long s_nnz = 0;
    for(int j=0; j<comm->Ns; j++){
      proc = comm->Nss[j];
      for(srow_id=0; srow_id<comm->S[proc]; srow_id++){
	if(gID == comm->LG[comm->SLID[proc][srow_id]]){
	  goto found_match;
	}
      }
    }

  found_match:
    for(int j=0; j<srow_id; j++){
      s_nnz += comm->SAp[proc][j];
    }
    for(int j=0; j<ncols; j++){
      gID = G_col_id[cols[j]];
      for(int k=0; k<comm->SAp[proc][srow_id]; k++){
	if(gID == comm->SGRId[proc][s_nnz+k]){
	  Lk[proc][s_nnz+k] += lk[send[i]*ncols +cols[j]];
	  break;
	}
      }
    }
  }/* end send loop */

  /* free remaining memory */
  free(cols);
  free(send);
  free(G_row_id);
  free(G_col_id);
}
