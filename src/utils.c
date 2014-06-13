/* HEADER */
#include "utils.h"
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include "mkl_cblas.h"
#include "mkl_lapack.h"

#include "enumerations.h"
#include "index_macros.h"
#include "get_dof_ids_on_elem.h"
#include "elem3d.h"
#include "homogen.h"
#include "allocation.h"
#include "incl.h"
#include "localizat.h"
#include "matice.h"
#include "resice.h"

#ifndef UTILS_DEBUG
#define UTILS_DEBUG 0
#endif

static const int periodic = 0;

int alloc_sprintf(char **str,
		  const char *format,
		  ...)
{
  int str_len = 0;
  va_list args;

  /* get length of string */
  va_start(args,format);
  str_len += vsnprintf(NULL,0,format,args);
  va_end(args);

  /* allocate and populate string */
  *str = PGFEM_calloc(str_len+1,sizeof(char));
  str_len = 0;
  va_start(args,format);
  str_len +=  vsprintf(*str,format,args);
  va_end(args);

  return str_len;
}

void pack_data(const void *src,
	       char *buffer,
	       size_t *pos,
	       const size_t n_el,
	       const size_t size)
{
  memcpy(buffer + *pos,src,n_el*size);
  *pos += n_el*size;
}

void unpack_data(const char *buffer,
		 void *dest,
		 size_t *pos,
		 const size_t n_el,
		 const size_t size)
{
  memcpy(dest,buffer + *pos,n_el*size);
  *pos += n_el*size;
}	

void mat2array(double *array,
	       const double **mat,
	       const unsigned int nrows,
	       const unsigned int ncols)
{
  for (int i=0; i<nrows; i++){
    memcpy(&array[idx_2_gen(i,0,nrows,ncols)],
	   mat[i],
	   ncols*sizeof(double));
  }
}


void array2mat(const double *array,
	       double **tensor,
	       const unsigned int I,
	       const unsigned int J)
{
  for (int i=0; i<I; i++){
    memcpy(tensor[i],
	   &array[idx_2_gen(i,0,I,J)],
	   J*sizeof(double));
  }

  if (UTILS_DEBUG){
    for (int i=0; i<I; i++){
      for(int j=0; j<J; j++){
	PGFEM_printf("%f ",tensor[i][j]);
      }
      PGFEM_printf("\n");
    }
    PGFEM_printf("\n\n");
  }/* DEBUG */

}

void tensor3_2array(double *array,
		    const double ***tensor,
		    const unsigned int I,
		    const unsigned int J,
		    const unsigned int K)
{
  for (int i=0; i<I; i++){
    for(int j=0; j<J; j++){
      memcpy(&array[idx_3_gen(i,j,0,I,J,K)],
	     tensor[i][j],
	     K*sizeof(double));
    }
  }

 if(UTILS_DEBUG){
    for(int i=0; i<I*J*K; i++){
      PGFEM_printf("%f ",array[i]);
    }
    PGFEM_printf("\n\n");
  }

}

void array2tensor3(const double *array,
		   double ***tensor,
		   const unsigned int I,
		   const unsigned int J,
		   const unsigned int K)
{
  for (int i=0; i<I; i++){
    for(int j=0; j<J; j++){
      memcpy(tensor[i][j],
	     &array[idx_3_gen(i,j,0,I,J,K)],
	     K*sizeof(double));
    }
  }

  if (UTILS_DEBUG){
    for (int i=0; i<I; i++){
      for(int j=0; j<J; j++){
	for(int k=0; k<K; k++){
	  PGFEM_printf("%f ",tensor[i][j][k]);
	}
	PGFEM_printf("\n");
      }
      PGFEM_printf("\n\n");
    }
  }/* DEBUG */

}

void tensor4_2array(double *array,
		    const double ****tensor,
		    const unsigned int I,
		    const unsigned int J,
		    const unsigned int K,
		    const unsigned int L)
{
  for (int i=0; i<I; i++){
    for(int j=0; j<J; j++){
      for(int k=0; k<K; k++){
	memcpy(&array[idx_4_gen(i,j,k,0,I,J,K,L)],
	       tensor[i][j][k],
	       L*sizeof(double));
      }
    }
  }

 if(UTILS_DEBUG){
    for(int i=0; i<I*J*K*L; i++){
      PGFEM_printf("%f ",array[i]);
    }
    PGFEM_printf("\n\n");
  }

}

void array2tensor4(const double *array,
		   double ****tensor,
		   const unsigned int I,
		   const unsigned int J,
		   const unsigned int K,
		   const unsigned int L)
{
  for (int i=0; i<I; i++){
    for(int j=0; j<J; j++){
      for(int k=0; k<K; k++){
	memcpy(tensor[i][j][k],
	       &array[idx_4_gen(i,j,k,0,I,J,K,L)],
	       L*sizeof(double));
      }
    }
  }

  if (UTILS_DEBUG){
    for (int i=0; i<I; i++){
      for(int j=0; j<J; j++){
	for(int k=0; k<K; k++){
	  for(int l=0; l<L; l++){
	    PGFEM_printf("%f ",tensor[i][j][k][l]);
	  }
	  PGFEM_printf("\n");
	}
	PGFEM_printf("\n\n");
      }
    }
  }/* DEBUG */

}

void shapeTensor2array(double *array,
		       const double ****ST,
		       const unsigned int nne)
{
  for (int a=0; a<nne; a++){
    for (int b=0; b<3; b++){
      for (int i=0; i<3; i++){
	for (int j=0; j<3; j++){
	  array[idx_4_gen(a,b,i,j,nne,3,3,3)] = ST[i][j][b][a];
	}
      }
    }
  }
}

void array2shapeTensor(const double *array,
		       double ****ST,
		       const unsigned int nne)
{
 for (int a=0; a<nne; a++){
    for (int b=0; b<3; b++){
      for (int i=0; i<3; i++){
	for (int j=0; j<3; j++){
	  ST[i][j][b][a] =  array[idx_4_gen(a,b,i,j,nne,3,3,3)];
	}
      }
    }
  }
}

double det2x2(const double *mat)
{
  return (mat[0]*mat[3] - mat[1]*mat[2]);
}

double det3x3(const double *mat)
{
  return ( ( mat[0]*(mat[4]*mat[8] - mat[5]*mat[7]) )
	   + mat[1]*(mat[5]*mat[6] - mat[3]*mat[8])
	   + mat[2]*(mat[3]*mat[7] - mat[4]*mat[6]) );
}

double det4x4(const double *mat)
{
  return (mat[1]*mat[11]*mat[14]*mat[4]
	  - mat[1]*mat[10]*mat[15]*mat[4]
	  - mat[11]*mat[13]*mat[2]*mat[4]
	  + mat[10]*mat[13]*mat[3]*mat[4]
	  - mat[0]*mat[11]*mat[14]*mat[5]
	  + mat[0]*mat[10]*mat[15]*mat[5]
	  + mat[11]*mat[12]*mat[2]*mat[5]
	  - mat[10]*mat[12]*mat[3]*mat[5]
	  - mat[1]*mat[11]*mat[12]*mat[6]
	  + mat[0]*mat[11]*mat[13]*mat[6]
	  + mat[1]*mat[10]*mat[12]*mat[7]
	  - mat[0]*mat[10]*mat[13]*mat[7]
	  - mat[15]*mat[2]*mat[5]*mat[8]
	  + mat[14]*mat[3]*mat[5]*mat[8]
	  + mat[1]*mat[15]*mat[6]*mat[8]
	  - mat[13]*mat[3]*mat[6]*mat[8]
	  - mat[1]*mat[14]*mat[7]*mat[8]
	  + mat[13]*mat[2]*mat[7]*mat[8]
	  + mat[15]*mat[2]*mat[4]*mat[9]
	  - mat[14]*mat[3]*mat[4]*mat[9]
	  - mat[0]*mat[15]*mat[6]*mat[9]
	  + mat[12]*mat[3]*mat[6]*mat[9]
	  + mat[0]*mat[14]*mat[7]*mat[9]
	  - mat[12]*mat[2]*mat[7]*mat[9]);
}

double getJacobian(const double *mat,
		   const int elem_id,
		   int *err)
{
  double J = det3x3(mat);
  *err  = 0;
  if (J <= 0.0){
    int err_rank = 0;
    PGFEM_Error_rank(&err_rank);
    PGFEM_printerr("[%d] ERROR: Negative Jacobian!"
	    " J = %.5e (element %d)\n",err_rank,J,elem_id);
    *err = 1;
  }
  return J;
}

int inv2x2(const double *mat,
	   double *mat_inv)
{
  double A = det2x2(mat);
  if(A != 0.0){
    mat_inv[0] = 1.0/A*mat[3];
    mat_inv[1] = -1.0/A*mat[1];
    mat_inv[2] = -1.0/A*mat[2];
    mat_inv[3] = 1.0/A*mat[0];
    return 0;
  } else return 1;
}

int inv3x3(const double *mat,
	   double *mat_inv)
{
  double A = det3x3(mat);

  if (fabs(A) >= 1.0e-15){
    /* source: Wikipedia */
    mat_inv[0] = 1/A*(mat[4]*mat[8] - mat[5]*mat[7]);
    mat_inv[1] = -1/A*(mat[1]*mat[8] - mat[2]*mat[7]);
    mat_inv[2] = 1/A*(mat[1]*mat[5] - mat[2]*mat[4]);

    mat_inv[3] = -1/A*(mat[3]*mat[8] - mat[5]*mat[6]);
    mat_inv[4] = 1/A*(mat[0]*mat[8] - mat[2]*mat[6]);
    mat_inv[5] = -1/A*(mat[0]*mat[5] - mat[2]*mat[3]);

    mat_inv[6] = 1/A*(mat[3]*mat[7] - mat[4]*mat[6]);
    mat_inv[7] = -1/A*(mat[0]*mat[7] - mat[1]*mat[6]);
    mat_inv[8] = 1/A*(mat[0]*mat[4] - mat[1]*mat[3]);
    return 0;
  } else return 1;
}

int inv4x4(const double *mat,
	   double *mat_inv)
{
  double A = det4x4(mat);

  if(A != 0.0){
  mat_inv[0] = (-mat[11]*mat[14]*mat[5]
		+ mat[10]*mat[15]*mat[5]
		+ mat[11]*mat[13]*mat[6]
		- mat[10]*mat[13]*mat[7]
		- mat[15]*mat[6]*mat[9]
		+ mat[14]*mat[7]*mat[9]);

  mat_inv[1] = (mat[1]*mat[11]*mat[14]
		- mat[1]*mat[10]*mat[15]
		- mat[11]*mat[13]*mat[2]
		+ mat[10]*mat[13]*mat[3]
		+ mat[15]*mat[2]*mat[9]
		- mat[14]*mat[3]*mat[9]);

  mat_inv[2] = (-mat[15]*mat[2]*mat[5]
		+ mat[14]*mat[3]*mat[5]
		+ mat[1]*mat[15]*mat[6]
		- mat[13]*mat[3]*mat[6]
		- mat[1]*mat[14]*mat[7]
		+ mat[13]*mat[2]*mat[7]);

  mat_inv[3] = (mat[11]*mat[2]*mat[5]
		- mat[10]*mat[3]*mat[5]
		- mat[1]*mat[11]*mat[6]
		+ mat[1]*mat[10]*mat[7]
		+ mat[3]*mat[6]*mat[9]
		- mat[2]*mat[7]*mat[9]);

  mat_inv[4] = (mat[11]*mat[14]*mat[4]
		- mat[10]*mat[15]*mat[4]
		- mat[11]*mat[12]*mat[6]
		+ mat[10]*mat[12]*mat[7]
		+ mat[15]*mat[6]*mat[8]
		- mat[14]*mat[7]*mat[8]);

  mat_inv[5] = (-mat[0]*mat[11]*mat[14]
		+ mat[0]*mat[10]*mat[15]
		+ mat[11]*mat[12]*mat[2]
		- mat[10]*mat[12]*mat[3]
		- mat[15]*mat[2]*mat[8]
		+ mat[14]*mat[3]*mat[8]);

  mat_inv[6] = (mat[15]*mat[2]*mat[4]
		- mat[14]*mat[3]*mat[4]
		- mat[0]*mat[15]*mat[6]
		+ mat[12]*mat[3]*mat[6]
		+ mat[0]*mat[14]*mat[7]
		- mat[12]*mat[2]*mat[7]);

  mat_inv[7] = (-mat[11]*mat[2]*mat[4]
		+ mat[10]*mat[3]*mat[4]
		+ mat[0]*mat[11]*mat[6]
		- mat[0]*mat[10]*mat[7]
		- mat[3]*mat[6]*mat[8]
		+ mat[2]*mat[7]*mat[8]);

  mat_inv[8] = (-mat[11]*mat[13]*mat[4]
		+ mat[11]*mat[12]*mat[5]
		- mat[15]*mat[5]*mat[8]
		+ mat[13]*mat[7]*mat[8]
		+ mat[15]*mat[4]*mat[9]
		- mat[12]*mat[7]*mat[9]);

  mat_inv[9] = (-mat[1]*mat[11]*mat[12]
		+ mat[0]*mat[11]*mat[13]
		+ mat[1]*mat[15]*mat[8]
		- mat[13]*mat[3]*mat[8]
		- mat[0]*mat[15]*mat[9]
		+ mat[12]*mat[3]*mat[9]);

  mat_inv[10] = (-mat[1]*mat[15]*mat[4]
		 + mat[13]*mat[3]*mat[4]
		 + mat[0]*mat[15]*mat[5]
		 - mat[12]*mat[3]*mat[5]
		 + mat[1]*mat[12]*mat[7]
		 - mat[0]*mat[13]*mat[7]);

  mat_inv[11] = (mat[1]*mat[11]*mat[4]
		 - mat[0]*mat[11]*mat[5]
		 + mat[3]*mat[5]*mat[8]
		 - mat[1]*mat[7]*mat[8]
		 - mat[3]*mat[4]*mat[9]
		 + mat[0]*mat[7]*mat[9]);

  mat_inv[12] = (mat[10]*mat[13]*mat[4]
		 - mat[10]*mat[12]*mat[5]
		 + mat[14]*mat[5]*mat[8]
		 - mat[13]*mat[6]*mat[8]
		 - mat[14]*mat[4]*mat[9]
		 + mat[12]*mat[6]*mat[9]);

  mat_inv[13] = (mat[1]*mat[10]*mat[12]
		 - mat[0]*mat[10]*mat[13]
		 - mat[1]*mat[14]*mat[8]
		 + mat[13]*mat[2]*mat[8]
		 + mat[0]*mat[14]*mat[9]
		 - mat[12]*mat[2]*mat[9]);

  mat_inv[14] = (mat[1]*mat[14]*mat[4]
		 - mat[13]*mat[2]*mat[4]
		 - mat[0]*mat[14]*mat[5]
		 + mat[12]*mat[2]*mat[5]
		 - mat[1]*mat[12]*mat[6]
		 + mat[0]*mat[13]*mat[6]);

  mat_inv[15] = (-mat[1]*mat[10]*mat[4]
		 + mat[0]*mat[10]*mat[5]
		 - mat[2]*mat[5]*mat[8]
		 + mat[1]*mat[6]*mat[8]
		 + mat[2]*mat[4]*mat[9]
		 - mat[0]*mat[6]*mat[9]);
  for(int i=0; i<16; i++) mat_inv[i] /= A;
  return 0;
  } else return 1;
}

int inverse(double const* A,
	    const int M,
	    double *A_I)
{
  if (M <= 0){return 1;}
  int info;
  int *iPerm;
  int lwork = M*M;
  double *work;

  info = 0;
  switch(M){
  case 1:
    if(A[0] == 0.0){
      info = 1;
      break;
    } else{
      A_I[0] = 1.0/A[0];
      break;
    }

  case 2: info = inv2x2(A,A_I); break;
  case 3: info = inv3x3(A,A_I); break;
  case 4: info = inv4x4(A,A_I); break;
  default:
    iPerm = aloc1i(M);
    work = aloc1(lwork);

    memcpy(A_I,A,M*M*sizeof(double));

    /* Factor into U matrix */
#ifdef ARCH_BGQ
    dgetrf(M,M,A_I,M,iPerm,&info);
#else
    dgetrf(&M,&M,A_I,&M,iPerm,&info);
#endif
    if(info<0){
      PGFEM_printerr("WARNING: illegal parameter given"
	      " to dgetrf at position %d.\n",info);
    } else if(info>0){
      PGFEM_printerr("WARNING: factor U is singular.\n");
    }

    /* Compute inverse using factored matrix */
#ifdef ARCH_BGQ
    dgetri(M,A_I,M,iPerm,work,lwork,&info);
#else
    dgetri(&M,A_I,&M,iPerm,work,&lwork,&info);
#endif
    if(info<0){
      PGFEM_printerr("WARNING: illegal parameter given"
	      " to dgetri at position %d.\n",info);
    } 
  
    free(iPerm);
    free(work);
    break;
  }/* switch M */

  if(info != 0){
    if(info > 0){
      PGFEM_printerr("ERROR: Matrix is singular, inverse not computed.\n");
    } else {
      PGFEM_printerr("ERROR: Error (%d) in inverse routine.\n",info);
    }
    PGFEM_Abort();
    abort();
  }

  return info;
}

void transpose(double *mat_t,
	       const double *mat,
	       const int mat_row,
	       const int mat_col)
{
  for(int i=0; i<mat_row; i++){
    for(int j=0; j<mat_col; j++){
      mat_t[idx_2_gen(j,i,mat_col,mat_row)] = 
	mat[idx_2_gen(i,j,mat_row,mat_col)];
    }
  }
}
	      
void symmetric_part(double *sym,
		    const double *mat,
		    const int dim)
{
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      sym[idx_2_gen(i,j,dim,dim)] = 0.5*(mat[idx_2_gen(i,j,dim,dim)]
					 + mat[idx_2_gen(j,i,dim,dim)]);
    }
  }
}

void print_coords(FILE *out,
		  const int nne,
		  const double *x,
		  const double *y,
		  const double *z)
{
  for(int i=0; i<nne; i++){
    PGFEM_fprintf(out,"%d: %22.15e %22.15e %22.15e\n",i,x[i],y[i],z[i]);
  }
  PGFEM_fprintf(out,"\n");
}

void print_array_d(FILE *out,
		   const double *array,
		   const int length,
		   const int nrow,
		   const int ncol)
{
  int n_block = length/(nrow*ncol);
  int count = 0;
  assert(n_block*nrow*ncol <= length);
  for (int i=0; i<n_block; i++){
    for (int j=0; j<nrow; j++){
      for (int k=0; k<ncol; k++){
	PGFEM_fprintf(out,"%12.5e ",array[count]);
	count++;
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n\n");
  }
}

void print_array_i(FILE *out,
		   const int *array,
		   const int length,
		   const int nrow,
		   const int ncol)
{
  int n_block = length/(nrow*ncol);
  int count = 0;
  assert(n_block*nrow*ncol <= length);
  for (int i=0; i<n_block; i++){
    for (int j=0; j<nrow; j++){
      for (int k=0; k<ncol; k++){
	PGFEM_fprintf(out,"%5d ",array[count]);
	count++;
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n\n");
  }
}

void print_array_l(FILE *out,
		   const long *array,
		   const int length,
		   const int nrow,
		   const int ncol)
{
  int n_block = length/(nrow*ncol);
  int count = 0;
  assert(n_block*nrow*ncol <= length);
  for (int i=0; i<n_block; i++){
    for (int j=0; j<nrow; j++){
      for (int k=0; k<ncol; k++){
	PGFEM_fprintf(out,"%5ld ",array[count]);
	count++;
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n\n");
  }
}

void print_material(FILE *out,
		    const MATERIAL *mat)
{
  PGFEM_fprintf(out,
	  "E_x,y,z:    %12.5e %12.5e %12.5e\n"
	  "G_yz,xz,xy: %12.5e %12.5e %12.5e\n"
	  "n_yz,xz,xy: %12.5e %12.5e %12.5e\n"
	  "a_x,y,z:    %12.5e %12.5e %12.5e\n"
	  "sig: %12.5e dev: %d vol: %d\n",
	  mat->Ex,mat->Ey,mat->Ez,
	  mat->Gyz,mat->Gxz,mat->Gxy,
	  mat->nyz,mat->nxz,mat->nxy,
	  mat->ax,mat->ay,mat->az,
	  mat->sig,mat->devPotFlag,mat->volPotFlag);
}

void update_elem_bub_dofs(const long ne,
			  ELEMENT *const elem)
{
  for(long i=0; i<ne; i++){
    cblas_daxpy(elem[i].n_bub*elem[i].n_bub_dofs,1.0,
		elem[i].d_bub_dofs,1,elem[i].bub_dofs,1);
    memset(elem[i].d_bub_dofs,0,
	   elem[i].n_bub*elem[i].n_bub_dofs*sizeof(double));
  }
}

void compute_disp_grad(const int nne,
		       const double *ST,
		       const double *disp,
		       double *grad,
		       const int node)
{
  /* Compute the gradient of a displacement field with the option to
     use a only a single node */
  memset(grad,0,9*sizeof(double));
  if(node >= 0 && node < nne){
    for(int b=0; b<3; b++){
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  grad[idx_2(i,j)] += (ST[idx_4_gen(node,b,i,j,nne,3,3,3)]
			       *disp[idx_2_gen(node,b,nne,3)]);
	}
      }
    }
  } else {
    for(int a=0; a<nne; a++){
      for(int b=0; b<3; b++){
	for(int i=0; i<3; i++){
	  for(int j=0; j<3; j++){
	    grad[idx_2(i,j)] += (ST[idx_4_gen(a,b,i,j,nne,3,3,3)]
				 *disp[idx_2_gen(a,b,nne,3)]);
	  }
	}
      }
    }
  }
}

long* list_boundary_el(const long ne,
		       const ELEMENT *elem,
		       const long nn,
		       const NODE *node,
		       const long myrank,
		       long *nbndel)
{
  long *nod;
  long *bnd;
  long *bndel;

  bnd = aloc1l(ne);
  *nbndel = 0;

  /* Determine the number of COMMUNICATION boundary elements and store
     bool for each element.  1 = boundary element */
  for(int i=0; i<ne; i++){
    nod = aloc1l(elem[i].toe);
    elemnodes(i,elem[i].toe,nod,elem);
    for(int j=0; j<elem[i].toe; j++){
      /*if(node[nod[j]].Gnn >= 0){*/
      if(node[nod[j]].Dom != myrank){
	bnd[i] = 1;
	(*nbndel)++;
	break;
      }
    }
    free(nod);
  }


  /* It is likely that nbndel << ne, therefore create lookup list for
     boundary elements */
  if(*nbndel == 0){
    PGFEM_printf("WARNING: no boundary elements on domain...seems odd\n");
    bndel = aloc1l(1);
  } else {
    bndel = aloc1l(*nbndel);
  }

  if(bndel == NULL){
    PGFEM_printf("Out of memory! (%s)\n",__func__);
    abort();
  }

  long count = 0;
  for(int i=0; i<ne; i++){
    if(bnd[i] == 1){
      bndel[count] = i;
      count++;
    }
  }

  free(bnd);

  if(count != *nbndel){
    PGFEM_printf("WARNING: number stored boundary elements != "
	   "number counted boundary elements!\n");
    PGFEM_printf("[%ld](%ld:%ld)\n", myrank,*nbndel,count);
  }

  return bndel;
}

long* times_print (FILE *in1,
		   const long nt,
		   const long n_p)
{
  long i,j,*print;
  double *help;
  
  if(nt == 0){
    print = aloc1l (1);
  } else {
    print = aloc1l (nt);
  }

  if(n_p == 0){
    help  = aloc1 (1);
  } else {
    help = aloc1 (n_p);
  }
  
  for (i=0;i<n_p;i++){
    fscanf (in1,"%lf",&help[i]);
  }
  
  for (i=0;i<n_p;i++){
    for (j=0;j<nt;j++){
      if (j == help[i]){
	print[j] = 1;
	continue;
      }
    }
  }
  
  dealoc1 (help);
  
  return (print);
}

long num_fib (long nmat,
	      long ne,
	      ELEMENT *elem)
/*
  funkce vraci pocet fiberu
*/
{
  long i,n,*a;
  
  a = (long*) PGFEM_calloc (nmat,sizeof(long));
  
  for (i=0;i<ne;i++) a[elem[i].mat[1]]=1;
  
  n=0;
  for (i=0;i<nmat;i++){
    if (a[i]==1)  n++;
  }
  
  free (a);
  return (n);
}

long num_matr (long nmat,
	       long ne,
	       ELEMENT *elem)
/*
  funkce vraci pocet matric
*/
{
  long i,n,*a;
  
  a = (long*) PGFEM_calloc (nmat,sizeof(long));
  
  for (i=0;i<ne;i++) a[elem[i].mat[0]] = 1;
  
  n=0;
  for (i=0;i<nmat;i++){
    if (a[i] == 1)  n++;
  }
  
  free (a);
  return (n);
}

long list (long ***a,
	   long ne,
	   long nmat,
	   long nc,
	   ELEMENT *elem)
/*
  function creates list of combinations fiber - matrix - volume fraction
*/
{
  long i,j,k,n;
  
  for (i=0;i<ne;i++){
    a[elem[i].mat[0]][elem[i].mat[1]][elem[i].hom[0]] = 1;
  }
  
  n=0;
  for (i=0;i<nmat;i++){
    for (j=0;j<nmat;j++){
      for (k=0;k<nc;k++){
	if (a[i][j][k] == 1) n++;
      }
    }
  }
  
  return (n);
}

double Tetra_V (const double *x,
		const double *y,
		const double *z)
/*
       
 */
{
  double V;
  
  V = 1./6.*(((x[1]-x[0])*(y[2]-y[0])*(z[3]-z[0]) +
	      (y[1]-y[0])*(z[2]-z[0])*(x[3]-x[0]) +
	      (z[1]-z[0])*(x[2]-x[0])*(y[3]-y[0]))-
	     ((z[1]-z[0])*(y[2]-y[0])*(x[3]-x[0]) +
	      (x[1]-x[0])*(z[2]-z[0])*(y[3]-y[0]) + 
	      (y[1]-y[0])*(x[2]-x[0])*(z[3]-z[0])));
  
  if (V < 0.0) V = -1.0*V;
  
  return (V);
}

double Tetra_qv_V (const long nne,
		   const long ndofn,
		   const double *x,
		   const double *y,
		   const double *z)
/*
       
 */
{
  long i,j,k,II,JJ,KK,ndofe;
  double *gk,*ge,*gz,*w,J,ksi,eta,zet,ai,aj,ak,**B_T,V;
  
  ndofe = nne*ndofn;
  
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1(5);
  B_T = aloc2(ndofe,6);
  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  V = 0.0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	
	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	V += ai*aj*ak*J;
	
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  dealoc2 (B_T,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  
  return (V);
}

double Hexa_V (const double *x,
	       const double *y,
	       const double *z)
{
  double *gk,*ge,*gz,*w;
  double J,ksi,eta,zet,ai,aj,ak,**B_T,V2;
  long i,j,k,II,JJ,KK,ndofe,nne;

  /* V = 0.0; */

  /* xx = aloc1(4); yy = aloc1(4); zz = aloc1(4); */

  /* /\* I *\/ */
  /* xx[0] = x[0];  yy[0] = y[0];  zz[0] = z[0]; */
  /* xx[1] = x[1];  yy[1] = y[1];  zz[1] = z[1]; */
  /* xx[2] = x[3];  yy[2] = y[3];  zz[2] = z[3]; */
  /* xx[3] = x[5];  yy[3] = y[5];  zz[3] = z[5]; */
  
  /* V +=  Tetra_V (xx,yy,zz); */

  /* /\* II *\/ */
  /* xx[0] = x[0];  yy[0] = y[0];  zz[0] = z[0]; */
  /* xx[1] = x[3];  yy[1] = y[3];  zz[1] = z[3]; */
  /* xx[2] = x[4];  yy[2] = y[4];  zz[2] = z[4]; */
  /* xx[3] = x[5];  yy[3] = y[5];  zz[3] = z[5]; */
  
  /* V +=  Tetra_V (xx,yy,zz); */
  
  /* /\* III *\/ */
  /* xx[0] = x[3];  yy[0] = y[3];  zz[0] = z[3]; */
  /* xx[1] = x[7];  yy[1] = y[7];  zz[1] = z[7]; */
  /* xx[2] = x[4];  yy[2] = y[4];  zz[2] = z[4]; */
  /* xx[3] = x[5];  yy[3] = y[5];  zz[3] = z[5]; */
  
  /* V +=  Tetra_V (xx,yy,zz); */
  
  /* /\* IV *\/ */
  /* xx[0] = x[1];  yy[0] = y[1];  zz[0] = z[1]; */
  /* xx[1] = x[2];  yy[1] = y[2];  zz[1] = z[2]; */
  /* xx[2] = x[3];  yy[2] = y[3];  zz[2] = z[3]; */
  /* xx[3] = x[7];  yy[3] = y[7];  zz[3] = z[7]; */
  
  /* V +=  Tetra_V (xx,yy,zz); */

  /* /\* V *\/ */
  /* xx[0] = x[1];  yy[0] = y[1];  zz[0] = z[1]; */
  /* xx[1] = x[2];  yy[1] = y[2];  zz[1] = z[2]; */
  /* xx[2] = x[5];  yy[2] = y[5];  zz[2] = z[5]; */
  /* xx[3] = x[7];  yy[3] = y[7];  zz[3] = z[7]; */
  
  /* V +=  Tetra_V (xx,yy,zz); */
  
  /* /\* VI *\/ */
  /* xx[0] = x[2];  yy[0] = y[2];  zz[0] = z[2]; */
  /* xx[1] = x[5];  yy[1] = y[5];  zz[1] = z[5]; */
  /* xx[2] = x[6];  yy[2] = y[6];  zz[2] = z[6]; */
  /* xx[3] = x[7];  yy[3] = y[7];  zz[3] = z[7]; */
  
  /* V +=  Tetra_V (xx,yy,zz); */
  
  /****************************************************/
  /****************************************************/
  /****************************************************/
  
  nne = 8; ndofe = nne*3;
  
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1(5);
  B_T = aloc2(ndofe,6);
  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  V2 = 0.0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){

	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
	}
	
	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	V2 += ai*aj*ak*J;
	
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  dealoc2 (B_T,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  
  /* dealoc1(xx); */
  /* dealoc1(yy); */
  /* dealoc1(zz); */
  
  return (V2);
}

double T_VOLUME (const long ne,
		 const long ndofn,
		 const ELEMENT *elem,
		 const NODE *node)
{
  long ii,nne,*nod;
  double PL=0.0,*x,*y,*z,volume;
  
  for (ii=0;ii<ne;ii++){
    
    nne = elem[ii].toe;
    
    nod = aloc1l (nne); x = aloc1 (nne); y = aloc1 (nne); z = aloc1 (nne);
    
    elemnodes (ii,nne,nod,elem);
    nodecoord_updated (nne,nod,node,x,y,z);
    /* Volume of element */
    if (nne == 4)   volume = Tetra_V (x,y,z);
    if (nne == 8)   volume = Hexa_V (x,y,z);
    if (nne == 10)  volume = Tetra_qv_V (nne,ndofn,x,y,z);
    PL += volume;
    
    dealoc1 (x); dealoc1 (y); dealoc1 (z); dealoc1l (nod);
  }
  
  return (PL);
}

double area (long nne,
	     double *x,
	     double *y)
/*
  Function returns area of element     
*/
{
  double a;
  
  a = 0.0;
  
  if (nne == 3){
    a = 0.5*(x[0]*y[1] + y[0]*x[2] + x[1]*y[2] - y[1]*x[2] - x[0]*y[2] - y[0]*x[1]);
  }
  if (nne == 4){
    a = 0.5*(x[0]*y[1] - x[1]*y[0] +
	     x[1]*y[2] - x[2]*y[1] +
	     x[2]*y[3] - x[3]*y[2] +
	     x[3]*y[0] - x[0]*y[3]);
  }
  return (a);
}

void def_elem (const long *cn,
	       const long ndofe,
	       const double *r,
	       const ELEMENT *elem,
	       const NODE *node,
	       double *r_e,
	       const SUPP sup,
	       const long TYPE)
{ 
  enum{UPDATED=0,TOTAL=1,OTHER=2};
  long i,j;

  for (i=0;i<ndofe;i++){
    j = cn[i];
    if (j == 0){
      r_e[i] = 0.0;
    } else if (j > 0){
      r_e[i] = r[j-1];
    } else { /* if (j < 0)  { */
      if (TYPE == UPDATED) r_e[i] = sup->defl_d[abs(j)-1];
      if (TYPE == TOTAL) r_e[i] = sup->defl[abs(j)-1];
      if (TYPE == OTHER) r_e[i] = 0.0;
    }
  }
}

void def_elem_total (const long *cn,
		     const long ndofe,
		     const double *r,
		     const double *d_r,
		     const ELEMENT *elem,
		     const NODE *node,
		     const SUPP sup,
		     double *r_e)
{
  for(int i=0; i< ndofe; i++){
    const int id = cn[i];
    const int aid = abs(id) - 1;

    if (id == 0){
      r_e[i] = 0.0;
    } else if (id > 0){
      r_e[i] = r[aid] + d_r[aid];
    } else {
      r_e[i] = sup->defl[aid] + sup->defl_d[aid];
    }
  }
}

void elemnodes (const long ii,
		const long nne,
		long *nod,
		const ELEMENT *elem)
/*
  returns nodes of actual element
       
  ii - index of actual element
  nne - number of nodes on element
  nod - array of nodes on element
       
  %%%%%%%%%%%%%%%% TESTED 7.12.99 %%%%%%%%%%%%%%%%%
*/
{
  long i;
  
  for (i=0;i<nne;i++){
    nod[i] = elem[ii].nod[i];
  }
}

/* void nodecoord (const long nne, */
/* 		const long *nod, */
/* 		const NODE *node, */
/* 		double *x, */
/* 		double *y, */
/* 		double *z) */
/* { */
/*   /\* Total Lagrangian *\/  */
/*   if (periodic == 1 || analysis == DISP){ */
/*     nodecoord_total(nne,nod,node,x,y,z); */
/*   } else { /\* Updated Lagrangian *\/ */
/*     nodecoord_updated(nne,nod,node,x,y,z); */
/*   } */
/* } */

void nodecoord_total (const long nne,
		      const long *nod,
		      const NODE *node,
		      double *x,
		      double *y,
		      double *z)
{
  for (int i=0;i<nne;i++){
    x[i] = node[nod[i]].x1_fd;
    y[i] = node[nod[i]].x2_fd;
    z[i] = node[nod[i]].x3_fd;
  }
}

void nodecoord_updated (const long nne,
			const long *nod,
			const NODE *node,
			double *x,
			double *y,
			double *z)
{
  for (int i=0;i<nne;i++){
    x[i] = node[nod[i]].x1;
    y[i] = node[nod[i]].x2;
    z[i] = node[nod[i]].x3;
  }
}


void list_el_prescribed_def (SUPP sup,
			     const NODE *node,
			     const ELEMENT *elem,
			     const BOUNDING_ELEMENT *b_elems,
			     const long ne,
			     const int n_be,
			     const long nn)
/*
       
 */
{
  long i,ii,j,k,nne,*nod,pom;
  long *key;
  val_key *lbnd_el;
  
  nod = aloc1l (10);
  
  if (sup->nde == 0)  sup->lepd = (long *) PGFEM_calloc (1,sizeof(long));
  else  sup->lepd = (long *) PGFEM_calloc (sup->nde,sizeof(long));

  if(n_be == 0){
    lbnd_el = PGFEM_calloc(1,sizeof(val_key));
    key = aloc1l(1);
  } else {
    lbnd_el = PGFEM_calloc(n_be,sizeof(val_key));
    key = aloc1l(n_be);
  }

  /* get the list of volume elements with prescribed deflection */
  ii = 0;
  for (i=0;i<ne;i++){
    nne = elem[i].toe;
    
    elemnodes (i,nne,nod,elem);
    
    pom = 1;
    for (j=0;j<sup->ndn;j++){
      for (k=0;k<nne;k++){
	if (sup->lnpd[j] == nod[k]){
	  pom = 0; break;
	}
      }
      if (pom == 0) break;
    }
    if (pom == 0) {
      sup->lepd[ii] = i;
      ii++;
    }
  }
  dealoc1l (nod);

  /* ADDED 1/7/2013 TESTED */
  /* there is now a filtered list of volume elements with prescribed
     deflections that is ordered by element index by construction. Get
     the volume element inidces that are on the bounding elements and
     sort them, keeping the associated bounding element index in
     tact. Then search the list of elements with prescribed
     deflections for a match and store. */
  for(i=0; i<n_be; i++){
    key[i] = -1;
    lbnd_el[i].key = i;
    lbnd_el[i].val = b_elems[i].vol_elem_id;
  }

  qsort(lbnd_el,n_be,sizeof(val_key),compare_val_w_key);

  ii = 0;
  sup->nd_be = 0;
  for(i=0; i<n_be; i++){
    for(j=ii; j<sup->nde; j++){
      /* found match. NOTE: ii is not incremented on match because
	 multiple bounding elements can bound the same volume
	 element */
      if(sup->lepd[j] == lbnd_el[i].val){
	key[i] = lbnd_el[i].key;
	sup->nd_be ++;
	break;
      } else if (sup->lepd[j] < lbnd_el[i].val){
	ii++;
      }
    }
  }

  /* the list boundiing elements with prescribed deflections is
     contained in 'key'. Sort the list and store the last sup->nd_be
     values. */
  if(sup->nd_be > 0){
    sup->lbepd = aloc1l(sup->nd_be);
    qsort(key,n_be,sizeof(long),compare_long);
    memcpy(sup->lbepd,(key + n_be - sup->nd_be),sup->nd_be*sizeof(long));
  } else {
    sup->lbepd = aloc1l(1);
  }

  free(lbnd_el);
  free(key);
}

void fun_eps (double *r_e,
	      long ndofe,
	      double **B_T,
	      double *eps)
/*
       
 */
{
  long i;
  double **r,**A;
  
  r = aloc2(1,ndofe); A = aloc2(1,6);
  
  for (i=0;i<ndofe;i++)  r[0][i] = r_e[i];
  nas_AB (r,B_T,A,1,ndofe,6);
  for (i=0;i<6;i++)  eps[i] = A[0][i];
  
  dealoc2(r,1); dealoc2(A,1);
}

void eps_element (long nne,
		  long ndofn,
		  double V,
		  double *r_e,
		  double *x,
		  double *y,
		  double *z,
		  double *EPSi)
/*
       
 */
{
  
  long i,j,k,II,JJ,KK,jj,ndofe;
  double *gk,*ge,*gz,*w,*eps,J,ksi,eta,zet,ai,aj,ak,**B_T;
  
  ndofe = nne*ndofn;
  
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1(5);
  eps = aloc1(6);
  B_T = aloc2(ndofe,6);

  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  nulld (EPSi,6);
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 4)  {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
	}
	
	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	/*  Matice B_BAR */
	B_BAR (B_T,nne,x,y,z);
	
	fun_eps (r_e,ndofe,B_T,eps);
	
	for (jj=0;jj<6;jj++)  EPSi[jj] += ai*aj*ak*eps[jj]*J/V;
	
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/

  dealoc2 (B_T,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (eps);  
}

void eps_e_in (long nne,
	       long ndofn,
	       double *r_e,
	       double *x,
	       double *y,
	       double *z,
	       double **EPSi)
/*
       
 */
{
  
  long i,j,k,II,JJ,KK,jj,ndofe,ip;
  double *gk,*ge,*gz,*w,*eps,J,ksi,eta,zet,ai,aj,ak,**B_T;
  
  ndofe = nne*ndofn;
  
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1(5);
  eps = aloc1(6);
  B_T = aloc2(ndofe,6);

  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  ip = 0;  nulld2 (EPSi,II*JJ*KK,6);
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 4)  {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
	}
	
	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	/*  Matice B_BAR */
	B_BAR (B_T,nne,x,y,z);
	
	fun_eps (r_e,ndofe,B_T,eps);
	
	for (jj=0;jj<6;jj++)  EPSi[ip][jj] = eps[jj];
	
	ip++;
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/

  dealoc2 (B_T,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (eps);  
}

void stress (long ne,
	     long ndofn,
	     NODE *node,
	     ELEMENT *elem,
	     MATGEOM matgeom,
	     HOMMAT *hommat,
	     double *r,
	     SIG *sig,
	     EPS *eps,
	     SUPP sup,
	     const int analysis)
/*
       
 */
{
  long ii,i,j,*nod,nne,ndofe,*cn;
  double *r_e,*x,*y,*z,*EPSi,**D,V;
  
  EPSi = aloc1(6); D = aloc2 (6,6);
  
  for (ii=0;ii<ne;ii++){
    
    nne = elem[ii].toe; 
    ndofe = nne*ndofn;
    
    r_e = aloc1(ndofe);
    x = aloc1(nne);
    y = aloc1(nne);
    z = aloc1(nne);
    nod = aloc1l (nne);
    cn = aloc1l (ndofe);

    
    /* vector of nodes on the element */
    elemnodes (ii,nne,nod,elem);
    /* nodal coordinates of the element */
    switch(analysis){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
    /* Id numbers */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    /* vector of nodal deformation on the element */
    def_elem (cn,ndofe,r,elem,node,r_e,sup,0);
    
    /* Volume of element */
    if (nne == 4)   V = Tetra_V (x,y,z);
    if (nne == 8)   V = Hexa_V (x,y,z);
    if (nne == 10)  V = Tetra_qv_V (nne,ndofn,x,y,z);
    
    /* epsilon on the element */	 
    eps_element (nne,ndofn,V,r_e,x,y,z,EPSi);
    
    /* material stiffnes matrix of the element */
    Stiffness_Matrix_3D (ii,0,elem,hommat,D,0);
    
    for (i=0;i<6;i++){
      eps[ii].el.o[i] = EPSi[i];
      for (j=0;j<6;j++){
	sig[ii].el.o[i] += D[i][j]*EPSi[j];
      }
    }
    dealoc1 (r_e);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1l (nod);
    dealoc1l (cn);
  }/*end of ii*/
  dealoc1 (EPSi); 
  dealoc2 (D,6);
}

void Mises_sig (long ne,
		SIG *sig,
		long TYPE)
/*
  Calculating of equvivalent Mises stresses vectors 
*/
{
  long i;
  double *Sij,Seq;
  
  Sij = aloc1(6);
  
  for (i=0;i<ne;i++){
    if (TYPE == 0){
      Sij[0] = sig[i].el.o[0] - (sig[i].el.o[0] + sig[i].el.o[1] + sig[i].el.o[2])/3.;
      Sij[1] = sig[i].el.o[1] - (sig[i].el.o[0] + sig[i].el.o[1] + sig[i].el.o[2])/3.;
      Sij[2] = sig[i].el.o[2] - (sig[i].el.o[0] + sig[i].el.o[1] + sig[i].el.o[2])/3.;
      Sij[3] = sig[i].el.o[3];
      Sij[4] = sig[i].el.o[4];
      Sij[5] = sig[i].el.o[5];
    }
    if (TYPE == 1){
      Sij[0] = sig[i].el.m[0] - (sig[i].el.m[0] + sig[i].el.m[1] + sig[i].el.m[2])/3.;
      Sij[1] = sig[i].el.m[1] - (sig[i].el.m[0] + sig[i].el.m[1] + sig[i].el.m[2])/3.;
      Sij[2] = sig[i].el.m[2] - (sig[i].el.m[0] + sig[i].el.m[1] + sig[i].el.m[2])/3.;
      Sij[3] = sig[i].el.m[3];
      Sij[4] = sig[i].el.m[4];
      Sij[5] = sig[i].el.m[5];
    }
    
    Seq = sqrt(3/2.*(Sij[0]*Sij[0] + Sij[1]*Sij[1] + Sij[2]*Sij[2] + 2*(Sij[3]*Sij[3] + Sij[4]*Sij[4] + Sij[5]*Sij[5])));
    
    if (TYPE == 0) sig[i].el.eq   = Seq;
    if (TYPE == 1) sig[i].el.eq_m = Seq;
  }
  dealoc1(Sij);
}

void Mises_eps (long ne,
		EPS *eps,
		long TYPE)
/*
       
 */
{
  long i;
  double *Eij,Eeq;
  
  Eij = aloc1(6);
  
  for (i=0;i<ne;i++){
    if (TYPE == 0){
      Eij[0] = eps[i].el.o[0] - (eps[i].el.o[0] + eps[i].el.o[1] + eps[i].el.o[2])/3.;
      Eij[1] = eps[i].el.o[1] - (eps[i].el.o[0] + eps[i].el.o[1] + eps[i].el.o[2])/3.;
      Eij[2] = eps[i].el.o[2] - (eps[i].el.o[0] + eps[i].el.o[1] + eps[i].el.o[2])/3.;
      Eij[3] = eps[i].el.o[3]/2.;
      Eij[4] = eps[i].el.o[4]/2.;
      Eij[5] = eps[i].el.o[5]/2.;
    }
    if (TYPE == 1){
      Eij[0] = eps[i].el.m[0] - (eps[i].el.m[0] + eps[i].el.m[1] + eps[i].el.m[2])/3.;
      Eij[1] = eps[i].el.m[1] - (eps[i].el.m[0] + eps[i].el.m[1] + eps[i].el.m[2])/3.;
      Eij[2] = eps[i].el.m[2] - (eps[i].el.m[0] + eps[i].el.m[1] + eps[i].el.m[2])/3.;
      Eij[3] = eps[i].el.m[3]/2.;
      Eij[4] = eps[i].el.m[4]/2.;
      Eij[5] = eps[i].el.m[5]/2.;
    }
    if (TYPE == 2){
      Eij[0] = eps[i].el.i[0] - (eps[i].el.i[0] + eps[i].el.i[1] + eps[i].el.i[2])/3.;
      Eij[1] = eps[i].el.i[1] - (eps[i].el.i[0] + eps[i].el.i[1] + eps[i].el.i[2])/3.;
      Eij[2] = eps[i].el.i[2] - (eps[i].el.i[0] + eps[i].el.i[1] + eps[i].el.i[2])/3.;
      Eij[3] = eps[i].el.i[3]/2.;
      Eij[4] = eps[i].el.i[4]/2.;
      Eij[5] = eps[i].el.i[5]/2.;
    }
    
    Eeq = sqrt (2/3.*(Eij[0]*Eij[0] + Eij[1]*Eij[1] + Eij[2]*Eij[2] + 2*(Eij[3]*Eij[3] + Eij[4]*Eij[4]+ Eij[5]*Eij[5])));
    
    if (TYPE == 0) eps[i].el.eq   = Eeq;
    if (TYPE == 1) eps[i].el.eq_m = Eeq;
    if (TYPE == 2) eps[i].el.eq_i = Eeq;
  }
  dealoc1 (Eij);
}

/********************************************/
/****** SMOOTHING OF STRESSES TO NODES ******/
/********************************************/

void stress_projector (long nne,
		       double *N,
		       double *P)
/*
       
 */
{
  long ii;
  double **A;
  
  A = aloc2 (nne,nne);
  
  nulld (P,nne);
  for (ii=0;ii<nne;ii++){
    P[ii] = N[ii];
  }
  dealoc2 (A,nne);
}

void str_prj_load (long ii,
		   long kk,
		   long nne,
		   long ndofn,
		   double *r_e,
		   double **D,
		   double *x,
		   double *y,
		   double *z,
		   SIG *sig_e,
		   double *f_e,
		   const int analysis)
/*
       
 */
{
  long i,j,jj,k,II,JJ,KK,ndofe,ip;
  double ksi,eta,zet,ai,aj,ak,J,**B_T,*gk,*ge,*gz,*w,*N,*eps,**EPSi,**sig,*P,*N_x,*N_y,*N_z;
  
  ndofe = ndofn*nne;
  
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);
  B_T = aloc2 (ndofe,6);
  N = aloc1 (nne);
  eps = aloc1 (6);
  EPSi = aloc2 (6,1);
  
  sig = aloc2 (6,1);
  P = aloc1 (nne);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  N_z = aloc1 (nne);
  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  nulld (f_e,nne);
  ip = 0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 4)  {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
	}

	shape_func (ksi,eta,zet,nne,N);
	
	switch(analysis){
	default:
	  J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
	  break;
	case ELASTIC:
	case TP_ELASTO_PLASTIC:
	  J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	  B_BAR (B_T,nne,x,y,z);
	  fun_eps (r_e,ndofe,B_T,eps);
	  for (jj=0;jj<6;jj++) EPSi[jj][0] = eps[jj];
	  nas_AB (D,EPSi,sig,6,6,1);
	  break;
	}
	
	stress_projector (nne,N,P);
	
	switch(analysis){
	default:
	  for (jj=0;jj<nne;jj++){
	    f_e[jj] += sig_e[ii].il[ip].o[kk]*ai*aj*ak*J*P[jj];
	  }
	  break;
	case ELASTIC:
	case TP_ELASTO_PLASTIC:
	  for (jj=0;jj<nne;jj++){
	    f_e[jj] += sig[kk][0]*ai*aj*ak*J*P[jj];
	  }
	  break;
	}
	
	ip++;
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  
  dealoc2 (B_T,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (P);
  dealoc2 (EPSi,6);
  dealoc2 (sig,6);
  dealoc1 (N);
  dealoc1 (eps);

  dealoc1 (N_x);
  dealoc1 (N_y);
  dealoc1 (N_z);
}

void str_solve (double *r,
		double *k,
		double *s,
		double *f,
		long *adr,
		long smo,
		long ne,
		long nn,
		long ndofn,
		NODE *node,
		ELEMENT *elem,
		HOMMAT *hommat,
		SIG *sig_e,
		SIG *sig_n,
		SUPP sup,
		const int analysis)
/*
       
 */
{
  
  long ii,i,j,*nod,nne,ndofe,*cn;
  double *r_e,*x,*y,*z,**D,*f_e;
  
  D = aloc2 (6,6);
  
  for (i=0;i<6;i++){
    nulld (f,smo);
    for (ii=0;ii<ne;ii++){
      
      nne = elem[ii].toe;  
      ndofe = nne*ndofn;
      
      r_e = aloc1 (ndofe);
      x = aloc1 (nne);
      y = aloc1 (nne);
      z = aloc1 (nne);
      nod = aloc1l (nne);
      f_e = aloc1 (nne);
      cn = aloc1l (ndofe);
      
      /* vector of nodes on the element */
      elemnodes (ii,nne,nod,elem);
      /* nodal coordinates of the element */
      switch(analysis){
      case DISP:
	nodecoord_total (nne,nod,node,x,y,z);
	break;
      default:
	nodecoord_updated (nne,nod,node,x,y,z);
	break;
      }

      /* Id numbers */
      get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
      
      if (analysis == ELASTIC || analysis == TP_ELASTO_PLASTIC) {
	/* material stiffnes matrix of the element */
	Stiffness_Matrix_3D (ii,0,elem,hommat,D,0);
	/* vector of nodal deformation on the element */
	def_elem (cn,ndofe,r,elem,node,r_e,sup,0);
      }
      
      str_prj_load (ii,i,nne,ndofn,r_e,D,x,y,z,sig_e,f_e,analysis);
      
      for (j=0;j<nne;j++){
	f[nod[j]] += f_e[j];
      }
      
      dealoc1l (nod);
      dealoc1 (r_e);
      dealoc1 (x);
      dealoc1 (y);
      dealoc1 (z);
      dealoc1 (f_e);
      dealoc1l (cn);
    }/*end ii < ne */
    
    nulld (s,smo);
    /*  solution of the system of equation  */
    resic_sky (k,s,f,adr,smo,3,analysis); /* Only back runing */
    
    for (j=0;j<nn;j++){
      sig_n[j].el.o[i] = s[j];
    }
  }/* end i < 6 */
  dealoc2 (D,6);
}

void str_elem_matrix (long kk,
		      long nne,
		      long ndofn,
		      double *x,
		      double *y,
		      double *z,
		      double *K)
/*
  For quadratic elements I need 11 points integration N*N
*/
{
  long i,ii,j,jj,k,II,JJ,KK,ndofe;
  double ksi,eta,zet,ai,aj,ak,J,**B_T,*gk,*ge,*gz,*w,*N,*P;
  
  ndofe = ndofn*nne;
  
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);
  B_T = aloc2 (ndofe,6);
  N = aloc1 (nne);
  P = aloc1 (nne);

  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  nulld (K,nne*nne);
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 4 || nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
	}

	shape_func (ksi,eta,zet,nne,N);
	
	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	stress_projector (nne,N,P);
	
	for (ii=0;ii<nne;ii++){
	  for (jj=0;jj<nne;jj++){
	    K[ii*nne+jj] += ai*aj*ak*J*P[ii]*N[jj];
	  }
	}
	
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  
  dealoc2 (B_T,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (N);
  dealoc1 (P);
}

void str_proj_matrix (long *adr,
		      long ne,
		      long ndofn,
		      ELEMENT *elem,
		      NODE *node,
		      HOMMAT *hommat,
		      double *k,
		      const int analysis)
/*
       
 */
{
  long i,j,nne,*cn,*nod;
  double *lk,*x,*y,*z,*s,*f;
  
  /***********************************************************************/
  /*  pro prvky s vice nez 10 uzly je treba pre delat nasledujici alokaci */
  /***********************************************************************/
  
  /* Google translation:
     for elements with more than 10 nodes is necessary to do the following
  */
  cn = aloc1l (10);
  nod = aloc1l (10);
  lk= aloc1 (10*10);
  x = aloc1 (10);
  y = aloc1 (10);
  z = aloc1 (10);
  
  s = aloc1 (10);
  f = aloc1 (10);
  
  for (i=0;i<ne;i++){
    
    nne = elem[i].toe;  
    elemnodes (i,nne,nod,elem);
    switch(analysis){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
    
    str_elem_matrix (i,nne,ndofn,x,y,z,lk);
    
    for (j=0;j<nne;j++) cn[j] = nod[j]+1;
    
    localizat (k,lk,adr,cn,nne);
  }
  dealoc1l (cn);
  dealoc1l (nod);
  dealoc1 (lk);
  dealoc1 (x);
  dealoc1(y);
  dealoc1 (z);
}

/********************************************************************************************/
/********************************************************************************************/

double eq_M_sig (long i,
		 long ip,
		 SIG *sig,
		 long TYPE)
/*
       
 */
{
  long k;
  double *Sij,*SS,*DS,S_eq;
  
  Sij = aloc1 (6);  SS = aloc1 (6); DS = aloc1 (6);
  
  for (k=0;k<6;k++){
    if (TYPE == 0) {SS[k] = sig[i].il[ip].m[k]; DS[k] = 0.0;}
    if (TYPE == 1) {SS[k] = sig[i].il[ip].m[k]; DS[k] = sig[i].d_il[ip].m[k];}
    if (TYPE == 2) {SS[k] = 0.0;                DS[k] = sig[i].d_il[ip].m[k];}
  }
  
  Sij[0] = (SS[0] + DS[0]) - ((SS[0] + DS[0]) + (SS[1] + DS[1]) + (SS[2] + DS[2]))/3.;
  Sij[1] = (SS[1] + DS[1]) - ((SS[0] + DS[0]) + (SS[1] + DS[1]) + (SS[2] + DS[2]))/3.;
  Sij[2] = (SS[2] + DS[2]) - ((SS[0] + DS[0]) + (SS[1] + DS[1]) + (SS[2] + DS[2]))/3.;
  Sij[3] = (SS[3] + DS[3]);
  Sij[4] = (SS[4] + DS[4]);
  Sij[5] = (SS[5] + DS[5]);
  
  S_eq = sqrt(3./2.*(Sij[0]*Sij[0] + Sij[1]*Sij[1] + Sij[2]*Sij[2] + 2*(Sij[3]*Sij[3] + Sij[4]*Sij[4] + Sij[5]*Sij[5])));
  
  dealoc1 (Sij); dealoc1 (SS); dealoc1 (DS);
  
  return (S_eq);
}

double eq_M_eps (long i,
		 long ip,
		 EPS *eps,
		 double **deps_i,
		 long TYPE)
/*
       
 */
{
  long k;
  double *Eij,*EE,*DE,E_eq;
  
  Eij = aloc1(6); EE = aloc1 (6); DE = aloc1 (6);
  
  for (k=0;k<6;k++){
    if (TYPE == 0) {EE[k] = eps[i].il[ip].i[k]; DE[k] = 0.0;}
    if (TYPE == 1) {EE[k] = eps[i].il[ip].i[k]; DE[k] = eps[i].d_il[ip].i[k] + deps_i[ip][k];}
    if (TYPE == 2) {EE[k] = 0.0;                DE[k] = eps[i].d_il[ip].i[k] + deps_i[ip][k];}
  }
  
  Eij[0] = (EE[0] + DE[0]) - ((EE[0] + DE[0]) + (EE[1] + DE[1]) + (EE[2] + DE[2]))/3.;
  Eij[1] = (EE[1] + DE[1]) - ((EE[0] + DE[0]) + (EE[1] + DE[1]) + (EE[2] + DE[2]))/3.;
  Eij[2] = (EE[2] + DE[2]) - ((EE[0] + DE[0]) + (EE[1] + DE[1]) + (EE[2] + DE[2]))/3.;
  Eij[3] = (EE[3] + DE[3])/2.;
  Eij[4] = (EE[4] + DE[4])/2.;
  Eij[5] = (EE[5] + DE[5])/2.;
  
  E_eq = sqrt(2/3.*(Eij[0]*Eij[0] + Eij[1]*Eij[1] + Eij[2]*Eij[2] + 2*(Eij[3]*Eij[3] + Eij[4]*Eij[4] + Eij[5]*Eij[5])));
  
  dealoc1 (Eij); dealoc1 (EE); dealoc1 (DE);
  
  return (E_eq);
}

void unequal_forces (long ii,
		     double *x,
		     double *y,
		     double *z,
		     long nne,
		     long ndofn,
		     ELEMENT *elem,
		     double **dsig,
		     double *fe)
/*
  f_int = BT*sig
*/
{
  long i,j,k,jj,kk,II,JJ,KK,ndofe,ip;
  double ksi,eta,zet,ai,aj,ak,J,**B_T;
  double *gk,*ge,*gz,*w;
  
  ndofe = ndofn*nne; nulld (fe,ndofe);
  
  gk = aloc1(5);  ge = aloc1(5);  gz = aloc1(5);  w = aloc1(5); B_T = aloc2(ndofe,6);
  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  ip = 0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){

	if (nne == 4)  {ksi = *(gk+k); eta = *(ge+k); zet = *(gz+k);  ai = *(w+k); aj = 1.0;    ak = 1.0;}
	if (nne == 10) {ksi = *(gk+k); eta = *(ge+k); zet = *(gz+k);  ai = *(w+k); aj = 1.0;    ak = 1.0;}
	if (nne == 8)  {ksi = *(gk+i); eta = *(gk+j); zet = *(gk+k);  ai = *(w+i); aj = *(w+j); ak = *(w+k);}

	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	/*  Matice B_BAR */
	B_BAR (B_T,nne,x,y,z);
	
	for (jj=0;jj<ndofe;jj++){
	  for (kk=0;kk<6;kk++){
	    fe[jj] += ai*aj*ak*J*B_T[jj][kk]*dsig[ip][kk];
	  }
	}
	ip++;
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  
  dealoc2 (B_T,ndofe);  dealoc1 (gk); dealoc1 (ge); dealoc1 (gz);  dealoc1 (w);
}

void aver_stress (long ii,
		  long nne,
		  long ndofn,
		  double *x,
		  double *y,
		  double *z,
		  SIG *sig,
		  EPS *eps)
/*
       
 */
{
  long i,j,k,II,JJ,KK,jj,ndofe,ip;
  double *gk,*ge,*gz,*w,J,ksi,eta,zet,ai,aj,ak,**B_T,V;
  
  ndofe = nne*ndofn;
  
  gk = aloc1(5);  ge = aloc1(5);  gz = aloc1(5);  w = aloc1(5);  B_T = aloc2(ndofe,6);
  
  /* Volume of element */
  if (nne == 4)   V = Tetra_V (x,y,z);
  if (nne == 8)   V = Hexa_V (x,y,z);
  if (nne == 10)  V = Tetra_qv_V (nne,ndofn,x,y,z);
  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  for (jj=0;jj<6;jj++){
    /* Stress */
    sig[ii].el.o[jj] = 0.0;
    sig[ii].el.m[jj] = 0.0;
    sig[ii].el.f[jj] = 0.0;
    /* Strain */
    eps[ii].el.o[jj] = 0.0;
    eps[ii].el.m[jj] = 0.0;
    eps[ii].el.f[jj] = 0.0;
    eps[ii].el.i[jj] = 0.0;
  }
  
  ip = 0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){

	if (nne == 4)  {ksi = *(gk+k); eta = *(ge+k); zet = *(gz+k);  ai = *(w+k); aj = 1.0;    ak = 1.0;}
	if (nne == 10) {ksi = *(gk+k); eta = *(ge+k); zet = *(gz+k);  ai = *(w+k); aj = 1.0;    ak = 1.0;}
	if (nne == 8)  {ksi = *(gk+i); eta = *(gk+j); zet = *(gk+k);  ai = *(w+i); aj = *(w+j); ak = *(w+k);}

	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	/*  Matice B_BAR */
	B_BAR (B_T,nne,x,y,z);
	
	for (jj=0;jj<6;jj++){
	  /* Stress */
	  sig[ii].el.o[jj] += ai*aj*ak*J/V*sig[ii].il[ip].o[jj];
	  sig[ii].el.m[jj] += ai*aj*ak*J/V*sig[ii].il[ip].m[jj];
	  sig[ii].el.f[jj] += ai*aj*ak*J/V*sig[ii].il[ip].f[jj];
	  /* Strain */
	  eps[ii].el.o[jj] += ai*aj*ak*J/V*eps[ii].il[ip].o[jj];
	  eps[ii].el.m[jj] += ai*aj*ak*J/V*eps[ii].il[ip].m[jj];
	  eps[ii].el.f[jj] += ai*aj*ak*J/V*eps[ii].il[ip].f[jj];
	  eps[ii].el.i[jj] += ai*aj*ak*J/V*eps[ii].il[ip].i[jj];
	}
	
	ip++;
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  dealoc2 (B_T,ndofe); dealoc1 (gk);  dealoc1 (ge);  dealoc1 (gz);  dealoc1 (w);
}

void check_equi (double *fu,
		 long ne,
		 long ndofd,
		 long ndofn,
		 ELEMENT *elem,
		 NODE *node,
		 MATGEOM matgeom,
		 SIG *sig,
		 const int analysis)
/*
       
 */
{
  long i,j,ii,II,JJ,*nod,nne,ndofe;
  double *x,*y,*z,**dsig,*fe;
  
  fe = aloc1 (30); nulld (fu,ndofd); 

  for (ii=0;ii<ne;ii++){
    
    nne = elem[ii].toe; ndofe = nne*ndofn;
    
    /* Integration */
    int_point (nne,&II);
    
    dsig = aloc2 (II,6); x = aloc1 (nne); y = aloc1 (nne); z = aloc1 (nne); nod = aloc1l (nne);
    
    /* vector of nodes on the element */
    elemnodes (ii,nne,nod,elem);
    /* nodal coordinates of the element */
    switch(analysis){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
    
    for (i=0;i<II;i++) for (j=0;j<6;j++)  dsig[i][j] = sig[ii].il[i].o[j] + sig[ii].d_il[i].o[j];
    
    /* Unequal_forces on the element */
    unequal_forces (ii,x,y,z,nne,ndofn,elem,dsig,fe);
    
    /* Localization */
    for (j=0;j<nne;j++){
      for (i=0;i<ndofn;i++){
	JJ = node[nod[j]].id[i]-1;
	if (JJ < 0)  continue;
	fu[JJ] += fe[j*ndofn+i];
      }/* end j */
    }/* end i */
    
    dealoc1 (x); dealoc1 (y); dealoc1 (z); dealoc1l (nod);
  }/* end elem */
  dealoc1 (fe);
}

double* Energy_functional (long ne,
			   long ndofn,
			   long ndofd,
			   ELEMENT *elem,
			   NODE *node,
			   SIG *sig,
			   EPS *eps,
			   MATGEOM matgeom,
			   double *f,
			   double *r,
			   const int analysis)
/*
       
 */
{
  long ii,i,j,k,jj,ip,nne,ndofe,II,JJ,KK,*nod;
  double *Sig,*Eps,*ENF,*gk,*ge,*gz,*w,J,ksi,eta,zet,ai,aj,ak,**B_T,*x,*y,*z,SE;
  
  
  ENF  = aloc1 (3);
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1(5);
  Sig = aloc1 (6);
  Eps = aloc1 (6);
  
  for (ii=0;ii<ne;ii++){
    nne = elem[ii].toe; ndofe = nne*ndofn;
    
    B_T = aloc2 (ndofe,6);
    
    /* Integration */
    integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
    
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
    nod = aloc1l (nne);
    
    /* vector of nodes on the element */
    elemnodes (ii,nne,nod,elem);
    /* nodal coordinates of the element */
    switch(analysis){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
    
    ip = 0;
    for (i=0;i<II;i++){
      for (j=0;j<JJ;j++){
	for (k=0;k<KK;k++){
	  
	  if (nne == 4)  {
	    ksi = *(gk+k);
	    eta = *(ge+k);
	    zet = *(gz+k);
	    ai = *(w+k);
	    aj = 1.0;
	    ak = 1.0;
	  }
	  if (nne == 10) {
	    ksi = *(gk+k);
	    eta = *(ge+k);
	    zet = *(gz+k);
	    ai = *(w+k);
	    aj = 1.0;
	    ak = 1.0;
	  }
	  if (nne == 8)  {
	    ksi = *(gk+i);
	    eta = *(gk+j);
	    zet = *(gk+k);
	    ai = *(w+i);
	    aj = *(w+j);
	    ak = *(w+k);
	  }
	  
	  J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	  
	  switch(analysis){
	  case FS_CRPL:
	  case FINITE_STRAIN:
	  case STABILIZED:
	    for (jj=0;jj<6;jj++){
	      Sig[jj] = sig[ii].il[ip].o[jj];
	      Eps[jj] = eps[ii].il[ip].o[jj];
	    }
	    break;
	  default:
	    for (jj=0;jj<6;jj++){
	      Sig[jj] = sig[ii].il[ip].o[jj] + sig[ii].d_il[ip].o[jj];
	      Eps[jj] = eps[ii].il[ip].o[jj] + eps[ii].d_il[ip].o[jj];
	    }
	    break;
	  }
	  
	  SE = ss (Sig,Eps,6);
	  
	  ENF[0] += ai*aj*ak*J*SE;
	  
	  ip++;
	}
      }
    }
    dealoc2 (B_T,ndofe);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1l (nod);
  }/* end elem */
  
  ENF[1] =  ss (f,r,ndofd);
  
  ENF[2] = ENF[0] - ENF[1];
  
  dealoc1 (Sig);
  dealoc1 (Eps);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
    
  return (ENF);
}

void tensor_9x9 (double **K,
		 double A[3][3][3][3],
		 long pom)
/*
       
 */
{
  long i,j,k,l,I,J;
  
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      if (i == 0 && j == 0) I = 0;
      if (i == 1 && j == 1) I = 1;
      if (i == 2 && j == 2) I = 2;

      if (i == 1 && j == 2) I = 3;
      if (i == 0 && j == 2) I = 4;
      if (i == 0 && j == 1) I = 5;

      if (i == 2 && j == 1) I = 6;
      if (i == 2 && j == 0) I = 7;
      if (i == 1 && j == 0) I = 8;
      
      for (k=0;k<3;k++){
	for (l=0;l<3;l++){
	  if (k == 0 && l == 0) J = 0;
	  if (k == 1 && l == 1) J = 1;
	  if (k == 2 && l == 2) J = 2;
 
	  if (k == 1 && l == 2) J = 3;
	  if (k == 0 && l == 2) J = 4;
	  if (k == 0 && l == 1) J = 5;

	  if (k == 2 && l == 1) J = 6;
	  if (k == 2 && l == 0) J = 7;
	  if (k == 1 && l == 0) J = 8; 
	  
	  if (pom == 0)  K[I][J] = A[i][j][k][l];
	  if (pom == 1)  A[i][j][k][l] = K[I][J];
	}
      }
    }
  }
}

double equivalent_Mises (long i,
			 SIG *sig)
/*
       
 */
{
  
  double *Sij,S_eq;
  
  Sij = aloc1(6);
  
  Sij[0] = (sig[i].el.o[0] - (sig[i].el.o[0] + sig[i].el.o[1]
			      + sig[i].el.o[2])/3.);

  Sij[1] = (sig[i].el.o[1] - (sig[i].el.o[0] + sig[i].el.o[1]
			      + sig[i].el.o[2])/3.);

  Sij[2] = (sig[i].el.o[2] - (sig[i].el.o[0] + sig[i].el.o[1]
			      + sig[i].el.o[2])/3.);

  Sij[3] = sig[i].el.o[3];
  Sij[4] = sig[i].el.o[4];
  Sij[5] = sig[i].el.o[5];
  
  S_eq = sqrt (3./2.*(Sij[0]*Sij[0] + Sij[1]*Sij[1] + Sij[2]*Sij[2]
		      + 2.*(Sij[3]*Sij[3] + Sij[4]*Sij[4]+ Sij[5]*Sij[5])));
  
  dealoc1 (Sij);
  
  return (S_eq);
}

double equivalent_M_eps (long i,
			 EPS *eps)
/*
       
 */
{
  
  double *Eij,E_eq;
  
  Eij = aloc1(6);
  
  Eij[0] = (eps[i].el.o[0] - (eps[i].el.o[0] + eps[i].el.o[1]
			      + eps[i].el.o[2])/3.);

  Eij[1] = (eps[i].el.o[1] - (eps[i].el.o[0] + eps[i].el.o[1]
			      + eps[i].el.o[2])/3.);

  Eij[2] = (eps[i].el.o[2] - (eps[i].el.o[0] + eps[i].el.o[1]
			      + eps[i].el.o[2])/3.);

  Eij[3] = eps[i].el.o[3]/2.;
  Eij[4] = eps[i].el.o[4]/2.;
  Eij[5] = eps[i].el.o[5]/2.;
  
  E_eq = sqrt (2./3.*(Eij[0]*Eij[0] + Eij[1]*Eij[1] + Eij[2]*Eij[2]
		      + 2.*(Eij[3]*Eij[3] + Eij[4]*Eij[4]+ Eij[5]*Eij[5])));
  
  dealoc1 (Eij);
  
  return (E_eq);
}

double equivalent_M_eps_pl (long i,
			    EPS *eps)
/*
       
 */
{
  
  double *Eij,E_eq;
  
  Eij = aloc1(6);
  
  Eij[0] = (eps[i].pl.o[0] - (eps[i].pl.o[0] + eps[i].pl.o[1]
			      + eps[i].pl.o[2])/3.);

  Eij[1] = (eps[i].pl.o[1] - (eps[i].pl.o[0] + eps[i].pl.o[1]
			      + eps[i].pl.o[2])/3.);

  Eij[2] = (eps[i].pl.o[2] - (eps[i].pl.o[0] + eps[i].pl.o[1]
			      + eps[i].pl.o[2])/3.);

  Eij[3] = eps[i].pl.o[3]/2.;
  Eij[4] = eps[i].pl.o[4]/2.;
  Eij[5] = eps[i].pl.o[5]/2.;
  
  E_eq = sqrt (2./3.*(Eij[0]*Eij[0] + Eij[1]*Eij[1] + Eij[2]*Eij[2]
		      + 2.*(Eij[3]*Eij[3] + Eij[4]*Eij[4]+ Eij[5]*Eij[5])));
  
  dealoc1 (Eij);
  
  return (E_eq);
}

void Mises (long ne,
	    SIG *sig,
	    EPS *eps,
	    const int analysis)
/*
  Calculating of equvivalent Mises stresses vectors 
*/
{
  long i;
  
  for (i=0;i<ne;i++){
    sig[i].el.eq = equivalent_Mises (i,sig);
    eps[i].el.eq = equivalent_M_eps (i,eps);
    if (analysis == FS_CRPL)
      eps[i].pl.eq[1] = equivalent_M_eps_pl (i,eps);
  }
}

void Logarithmic_strain (double **F,
			 double **EL)
/*
  v - eigen vectors
  eig - eigen values :: eig = lam^2
  lam - stretches in the principal directions :: lam = sqrt (eig)
*/
{
  double **C,**v,*e,*eig,*lam;
  int i,j,k;	
  
  eig = aloc1 (4);
  e = aloc1 (4);
  v = aloc2 (4,4);
  C = aloc2 (3,3);
  lam = aloc1 (3);
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	C[i][j] += F[k][i]*F[k][j];
  
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++){
      EL[i-1][j-1] = 0.0;
      if (fabs(C[i-1][j-1]) < 1.e-15)
	v[i][j] = 0.0;
      else
	v[i][j] = C[i-1][j-1];
    }
  
  /* Numerical recepises subroutines */
  tred2 (v,3,eig,e);
  tqli (eig,e,3,v);
  
  /* The unique eigenvalues are the squares of the stretches in the
     principal directions */
  lam[0] = sqrt(eig[1]);
  lam[1] = sqrt(eig[2]);
  lam[2] = sqrt(eig[3]);
 
  for(i=1;i<=3;i++){
    for(j=1;j<=3;j++){
      for(k=1;k<=3;k++){
	EL[j-1][k-1] += log(lam[i-1])*v[j][i]*v[k][i];
      }
    }
  }
  
  dealoc1(eig);
  dealoc1(e);
  dealoc2(v,4);
  dealoc2(C,3);
  dealoc1 (lam);
}

/*************************************************************************
 *               BEGIN FUNCTIONS UNIQUE TO PARALLEL CODE                 *
 *************************************************************************/

void LToG (const double *f,
	   double *Gf,
	   const int myrank,
	   const int nproc,
	   const long ndofd,
	   const long *DomDof,
	   const long GDof,
	   const COMMUN comm,
	   const MPI_Comm mpi_comm)
/*
       
 */
{
  long i,j,KK;
  double **send,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;
  
  /* for (i=0;i<DomDof[myrank];i++) */
  /*   Gf[i] = 0.0; */
  nulld(Gf,DomDof[myrank]);
  
  /* Allocate recieve */
  send = (double**) PGFEM_calloc (nproc,sizeof(double*));
  recieve = (double**) PGFEM_calloc (nproc,sizeof(double*));
  for (i=0;i<nproc;i++) {
    if (comm->S[i] == 0) KK = 1; else KK = comm->S[i];
    send[i] = (double*) PGFEM_calloc (KK,sizeof(double));
    if (comm->R[i] == 0) KK = 1; else KK = comm->R[i];
    recieve[i] = (double*) PGFEM_calloc (KK,sizeof(double));
  }

  /* Allocate status and request fields */
  if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
  sta_s = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
  req_s = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));

  if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
  sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
  req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));

  if(UTILS_DEBUG){/* debug the send-recieve */
    for(int i=0; i<nproc; i++){
      if(myrank == i){
	for(int j=0; j<nproc; j++){
	  PGFEM_printf("[%d]: sending %ld to %d\n",myrank,comm->S[j],j);
	  PGFEM_printf("[%d]: recieving %ld from %d\n",myrank,comm->R[j],j);
	}
      }
      MPI_Barrier(mpi_comm);
    }
    MPI_Barrier(mpi_comm);
  }


  /****************/
  /* Recieve data */
  /****************/
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    MPI_Irecv (recieve[KK],comm->R[KK],MPI_DOUBLE,KK,
	       MPI_ANY_TAG,mpi_comm,&req_r[i]);
  }

  /*************/
  /* Send data */
  /*************/
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    for (j=0;j<comm->S[KK];j++)
      send[KK][j] = f[comm->SLID[KK][j]];

    MPI_Isend (send[KK],comm->S[KK],MPI_DOUBLE,KK,
		myrank,mpi_comm,&req_s[i]);
  }

  pgfem_comm_get_owned_global_dof_values(comm,f,Gf);

  /* Wait to complete the comunicatins */
  MPI_Waitall (comm->Nr,req_r,sta_r);
  
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    for (j=0;j<comm->R[KK];j++)
      Gf[comm->RGID[KK][j] - GDof] += recieve[KK][j];
  }
  
  MPI_Waitall (comm->Ns,req_s,sta_s);
  /* Deallocate send and recieve */
  for (i=0;i<nproc;i++) {
    free(send[i]);
    free (recieve[i]);
  }
  free(send);
  free (recieve);

  free(sta_s);
  free(sta_r);	
  free(req_s);
  free(req_r);   
}

void GToL (const double *Gr,
	   double *r,
	   const int myrank,
	   const int nproc,
	   const long ndofd,
	   const long *DomDof,
	   const long GDof,
	   const COMMUN comm,
	   const MPI_Comm mpi_comm)
/*
       
 */
{
  long i,j,KK;
  double **send,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  //  for (i=0;i<ndofd;i++) r[i] = 0.0;
  nulld(r,ndofd);

  /* Allocate recieve */
  send = (double**) PGFEM_calloc (nproc,sizeof(double*));
  recieve = (double**) PGFEM_calloc (nproc,sizeof(double*));
  for (i=0;i<nproc;i++) {
    if (comm->R[i] == 0) KK = 1; else KK = comm->R[i];
    send[i] = (double*) PGFEM_calloc (KK,sizeof(double));
    if (comm->S[i] == 0) KK = 1; else KK = comm->S[i];
    recieve[i] = (double*) PGFEM_calloc (KK,sizeof(double));
  }
  
  /* Allocate status and request fields */
  if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
  sta_s = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
  req_s = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));

  if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
  sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
  req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));


  /****************/
  /* Recieve data */
  /****************/
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    MPI_Irecv (recieve[KK],comm->S[KK],MPI_DOUBLE,KK,
	       MPI_ANY_TAG,mpi_comm,&req_r[i]);
  }

  /*************/
  /* Send data */
  /*************/
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    for (j=0;j<comm->R[KK];j++)
      send[KK][j] = Gr[comm->RGID[KK][j] - GDof];

    MPI_Isend (send[KK],comm->R[KK],MPI_DOUBLE,KK,
	       myrank,mpi_comm,&req_s[i]);
  }
  
  pgfem_comm_get_local_dof_values_from_global(comm,Gr,r);

  /* Wait to complete the comunicatins */
  MPI_Waitall (comm->Ns,req_r,sta_r);
  
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    for (j=0;j<comm->S[KK];j++)
      r[comm->SLID[KK][j]] = recieve[KK][j];
  }
  
  MPI_Waitall (comm->Nr,req_s,sta_s);
  /* Deallocate send and recieve */
  for (i=0;i<nproc;i++) {
    free(send[i]);
    free (recieve[i]);
  }

 free(send);
 free (recieve);

 free(sta_s);
 free(sta_r);	
 free(req_s);
 free(req_r);   

}

MPI_Comm* CreateGraph (int nproc,
		       int myrank,
		       long nn,
		       NODE *node,
		       MPI_Comm mpi_comm)
/*
       
 */
{
  int *BN,*displ;
  long i,j,k,NBn=0,*hu1,TBn=0,*GNn,pom,Dom,II,*CDom;
  MPI_Comm *GrComm = NULL;
  
  for (i=0;i<nn;i++) if (node[i].Gnn >= 0) NBn++;
  
  hu1 = aloc1l (NBn);
  
  j = 0;
  for (i=0;i<nn;i++) if (node[i].Gnn >= 0) {hu1[j] = node[i].Gnn; j++;}
  
  /* Sort boundary nodes */
  for (j=0;j<NBn;j++){
    for (i=0;i<NBn-1;i++){
      while (hu1[i] > hu1[i+1]){
	pom = hu1[i]; hu1[i] = hu1[i+1]; hu1[i+1] = pom;
      }
    }
  }

  /* Gather number of domain boundary nodes */
  BN = aloc1i (nproc);
  MPI_Allgather (&NBn,1,MPI_LONG,BN,1,MPI_LONG,mpi_comm);

  for (i=0;i<nproc;i++) TBn += BN[i];
  displ = aloc1i (nproc);
  for (i=1;i<nproc;i++) displ[i] = displ[i-1] + BN[i-1];
  
  /* Gather nodes on boundaries between domains */
  GNn = aloc1l (TBn);
  MPI_Allgatherv (hu1,NBn,MPI_LONG,GNn,BN,displ,MPI_LONG,mpi_comm);
  
  for (i=0;i<nproc;i++){
    II = 0;
    if (i == myrank) continue;
    for (j=0;j<NBn;j++){
      for (k=0;k<BN[i];k++){
	if (hu1[j] == GNn[displ[i]+k]) {Dom++; II = 1; break;}
      }/* k < BN[i] */
      if (II == 1) break;
    }/* j < NBn */
  }/* i < nproc */

  CDom = aloc1l (Dom);
  
  Dom = 0;
  for (i=0;i<nproc;i++){
    II = 0;
    if (i == myrank) continue;
    for (j=0;j<NBn;j++){
      for (k=0;k<BN[i];k++){
	if (hu1[j] == GNn[displ[i]+k]) {CDom[Dom] = i; Dom++; II = 1; break;}
      }/* k < BN[i] */
      if (II == 1) break;
    }/* j < NBn */
  }/* i < nproc */
  
  /* PRINT */
  PGFEM_printf ("DOMAIN [%ld] is connected to ::",i);
  for (j=0;j<Dom;j++) PGFEM_printf ("  %ld",CDom[j]);
  PGFEM_printf ("\n");
  
  dealoc1l (hu1); dealoc1i (BN); dealoc1i (displ); dealoc1l (GNn); dealoc1l (CDom);

  return (GrComm);
}

void pause_time(int t)
{
  clock_t end;
  end = clock()+t*CLOCKS_PER_SEC;
  while(clock()<end){}
}

long* change_length(long *orig,
		    const long old_len,
		    const long new_len)
{
  long *pom,i;
  pom = aloc1l(new_len);
  if(new_len >= old_len) {for(i=0;i<old_len;i++) pom[i] = orig[i];}
  else {for(i=0;i<new_len;i++) pom[i] = orig[i];}
  free(orig);
  return pom;
}

void null_quit(void *array,
	       int error)
{
  if(array == NULL){
    PGFEM_printf("\nMemory full.\n");
    fflush(stdout);
    PGFEM_printf("Exit error %d\n",error);
    PGFEM_Abort();
  }
}

/*******************************************************************************************/
/*****************************  ARC-LENGTH PROCEDURES  *************************************/
/*******************************************************************************************/

long diag_K (double *k,
	     long *adr,
	     long ndofd)
/*
  P_D = +1 : Positive definite Matrix
  P_D = -1 : Negative definite Matrix
*/
{
  long i,P_D;
  P_D = 1;
  
  for (i=0;i<ndofd;i++){
    if (k[adr[i]] <= 0.0) P_D = -1;
  }
  
  return (P_D);
}

double det_K (double *k,
	      long *adr,
	      long ndofd)
/*
       
 */
{
  long i;
  double DET;
  
  DET = 1.0;
  /* for (i=0;i<ndofd;i++) DET *= k[adr[i]]; */
  for (i=0;i<ndofd;i++) DET *= k[adr[i]]/fabs(k[adr[i]]);
  
  return (DET);
}

double new_arc_length (long iter,
		       long iter_des,
		       double dAL,
		       double dAL0)
/*
       
 */
{
  double I1,I2,NEW;
  
  I1 = (double) iter_des;
  I2 = (double) iter;
  
  NEW = dAL*sqrt (I1/I2);
  if (NEW > dAL0) NEW = dAL0;
  
  return (NEW);
}


/***********************************************************/
/************  NONSYMMETRIC SPARSE SOLVER  *****************/
/***********************************************************/

long* sparse_ApAi (long ne,
		   long ndofd,
		   long ndofn,
		   ELEMENT *elem,
		   NODE *node,
		   long *Ap)
/*
  Sparse nonsymmetric column storage format Ap
*/
{
  long i,j,k,II,JJ,*cn,nne,ndofe,*nod,**AA,*ap,*Ai;
  
  /* Alocation has to be changed for elements with more than 10 nodes
     and 30 degreess of freedom */
  cn = aloc1l (30); nod = aloc1l (10); ap = aloc1l (ndofd);
  
  for (i=0;i<ne;i++){/* Number of contributions from elements to Aij */
    
    nne = elem[i].toe;  
    ndofe = ndofn*nne;
    elemnodes (i,nne,nod,elem);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    for (j=0;j<ndofe;j++){
      II = cn[j]-1;
      if (II < 0)  continue;
      for (k=0;k<ndofe;k++){
	JJ = cn[k]-1;
	if (JJ < 0)  continue;
	Ap[II]++;
      }
    }
  }/* end i < ne */
  
  AA = (long**) PGFEM_calloc (ndofd,sizeof(long*));
  for (i=0;i<ndofd;i++) AA[i]= (long*) PGFEM_calloc (Ap[i],sizeof(long));
  if (AA == NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    fflush(stdout);
    abort ();
  }
  
  for (i=0;i<ne;i++){/* List of ID number for Aij */
    
    nne = elem[i].toe;  
    ndofe = ndofn*nne;
    elemnodes (i,nne,nod,elem);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    for (j=0;j<ndofe;j++){
      II = cn[j]-1;
      if (II < 0)  continue;
      for (k=0;k<ndofe;k++){
	JJ = cn[k]-1;
	if (JJ < 0)  continue;
	AA[II][ap[II]] = JJ;
	ap[II]++;
      }
    }
  }/* end i < ne */
  
  for (k=0;k<ndofd;k++){/* Sort list of IDs for Aij */
    Ap[k] = 0;
    for (j=0;j<ap[k];j++){
      for (i=0;i<ap[k]-1;i++){
	while (AA[k][i] > AA[k][i+1]){
	  II = AA[k][i];
	  AA[k][i] = AA[k][i+1];
	  AA[k][i+1] = II;
	}
      }
    }
  }
  
  Ap[0] = 0;
  for (k=0;k<ndofd;k++){ /* Number of non-zeros in rows */
    for (j=0;j<ap[k]-1;j++){
      if (AA[k][j] < AA[k][j+1]) Ap[k+1]++;
    }
    Ap[k+1] += Ap[k] + 1;
  }
  
  Ai = aloc1l (Ap[ndofd]);
  
  for (k=0;k<ndofd;k++){/* Row indexes */
    i = Ap[k];
    for (j=0;j<ap[k]-1;j++){
      if (AA[k][j] < AA[k][j+1]) {Ai[i] = AA[k][j]; i++;}
    }
    Ai[i] = AA[k][j];
  }
  
  for (i=0;i<ndofd;i++) free (AA[i]); free (AA);
  dealoc1l (cn);  dealoc1l (nod); dealoc1l (ap);
  
  return (Ai);
}
