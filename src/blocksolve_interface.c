/* HEADER */

#include "blocksolve_interface.h"
#include "PGFEM_io.h"
#include "allocation.h"

void null_BSspmat (BSspmat *K)
/*
       
 */
{
  long i,j;
  
  for (i=0;i<K->num_rows;i++) {
    for (j=0;j<K->rows[i]->length;j++) {
      K->rows[i]->nz[j] = 0.0;
    }
  }
}

void write_mat_matlab(char *str,
		      BSspmat *A,
		      BSprocinfo  *procinfo)
{
  FILE	*fp;
  int	i, j, row, count = 1, nnz;
  
  /* try to open the file */
  if ((fp = fopen(str,"w")) == NULL) {
    PGFEM_printerr("Cannot open %s\n",str);
    exit(-1);
  }
  
  nnz = 0;
  for (i=0;i<A->num_rows;i++) {
    for (j=0;j<A->rows[i]->length;j++) {
      nnz++;
    }
  }
  
  /* print out the header */
  PGFEM_fprintf(fp,"n = %d;number_nz=%d;t=zeros(number_nz,3);\n",A->num_rows,nnz);
  
  /* print out the triplets */
  for (i=0;i<A->num_rows;i++) {
    A->map->flocal2global(1,&i,&row,procinfo,A->map);
    for (j=0;j<A->rows[i]->length;j++) {
      PGFEM_fprintf(fp,"t(%d,1:3) = [%d %d %4.16e];\n",
	      count,row+1,A->rows[i]->col[j]+1,A->rows[i]->nz[j]);
      count++;
    }
  }
  PGFEM_fprintf(fp,"a = sparse(t(:,1),t(:,2),t(:,3));\n");
  
  fclose(fp);
}

BSspmat *BSalloc_A (int start_num,
		    int n,
		    int *rp,
		    int *cval,
		    BSprocinfo *procinfo)
/* 
   This rutine is te same as BSeasy_A, but it onlu allocate the matrix 
*/
{
  BSspmat *A;
  int	i, j;
  int	*i_ptr;
  int	gnum;
  
  /* set up the structure and mapping for the sparse matrix */
  /* allocate the pointer to A */
  A = (BSspmat *) PGFEM_calloc (1,sizeof(BSspmat));
  
  /* set the number of local rows */
  A->num_rows = n;
  
  /* set the number of global rows */
  A->global_num_rows = n;
  GISUM(&(A->global_num_rows),1,&i,procinfo->procset);
  
  /* allocate the array of rows, and the space in each row */
  /* allow for a max length in each row, but set the current length to 0 */
  A->rows = (BSsprow **) PGFEM_calloc (A->num_rows,sizeof(BSsprow *));
  for (i=0;i<A->num_rows;i++) {
    A->rows[i] = (BSsprow *) PGFEM_calloc (1,sizeof(BSsprow));
    A->rows[i]->length = rp[i+1] - rp[i];
    A->rows[i]->col = &(cval[rp[i]]);
    A->rows[i]->nz = (double *) PGFEM_calloc (A->rows[i]->length,sizeof(double));
    gnum = start_num + i;
    A->rows[i]->diag_ind = -1;
    for (j=0;j<A->rows[i]->length;j++) {
      if (A->rows[i]->col[j] == gnum) {
	A->rows[i]->diag_ind = j;
	break;
      }
    }
  }
  
  /* allocate a pointer to a mapping structure */
  A->map = (BSmapping *) PGFEM_calloc (1,sizeof(BSmapping));
  
  /* set up the local to global mapping */
  /* all we need for this is the beginning number (offset) of */
  /* the local rows in the global numbering (see BSloc2glob) */
  A->map->vlocal2global = (void *) PGFEM_calloc (1,sizeof(int));
  i_ptr = (int *) A->map->vlocal2global; /* pointer to mapping data */
  *(i_ptr) = start_num;
  A->map->flocal2global = BSloc2glob; /* the mapping function */
  A->map->free_l2g = BSfreel2g; /* the routine to free the mapping */
  
  /* set up the global to local mapping */
  /* all we need for this is the beginning number (offset) of */
  /* the local rows in the global numbering (see BSglob2loc) */
  A->map->vglobal2local = (void *) PGFEM_calloc (1,sizeof(int));
  i_ptr = (int *) A->map->vglobal2local;	/* pointer to mapping data */
  *(i_ptr) = start_num;
  A->map->fglobal2local = BSglob2loc; /* the mapping function */
  A->map->free_g2l = BSfreeg2l;	/* the routine to free the mapping */
  
  /* set up the global to processor number mapping */
  /* we call the routine BSmake_off_map to create the mapping data */
  /* the local rows in the global numbering (see BSglob2proc) */
  A->map->vglobal2proc = (void *) BSmake_off_map(start_num,procinfo,A->global_num_rows); CHKERRN(0);
  A->map->fglobal2proc = BSglob2proc;	/* the mapping function */
  A->map->free_g2p = BSfree_off_map;	/* the routine to free the mapping */
  
  /* check for errors in A */
  if (procinfo->error_check) {
    BSrow_err_check(A,procinfo); CHKERRN(0);
  }
  return(A);
}
