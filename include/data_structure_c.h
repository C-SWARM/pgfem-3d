#ifndef _datastructure_c_H_
#define _datastructure_c_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mkl_cblas.h"

#define Define_Matrix(T)                                                \
typedef struct Matrix_##T                                               \
{                                                                       \
  long m_row, m_col;                                                    \
  long sizeof_T;                                                        \
  T *m_pdata;                                                           \
  T *temp;                                                              \
} Matrix_##T                                                            \


#define Matrix(T) Matrix_##T
#define Mat_v(p, m, n) p.m_pdata[(m-1)*p.m_col+(n-1)]

#define Matrix_construct(T, p) do {                                     \
  p.m_pdata = NULL;                                                     \
  p.temp    = NULL;                                                     \
  p.m_row = 0;                                                          \
  p.m_col = 0;                                                          \
  p.sizeof_T = sizeof(T);                                               \
} while(0) 

#define Matrix_construct_redim(T, p, m, n) do {                         \
  p.temp    = NULL;                                                     \
  p.m_row = m;                                                          \
  p.m_col = n;                                                          \
  p.sizeof_T = sizeof(T);                                               \
  p.m_pdata = malloc(p.sizeof_T*m*n);                                   \
} while(0)

#define Matrix_redim(p, m, n) do {                                      \
  Matrix_cleanup(p);                                                    \
  p.m_pdata = malloc(p.sizeof_T*m*n);                                   \
  p.temp    = NULL;                                                     \
  p.m_row = m;                                                          \
  p.m_col = n;                                                          \
} while(0)                                     

#define Matrix_init(p, value) do {                                      \
  long a, b;                                                            \
  for(a = 0; a < p.m_row*p.m_col; a++){                                 \
      p.m_pdata[a] =  value;                                            \
  }                                                                     \
} while(0) 

#define Matrix_construct_init(T, p, m, n, value) do {                   \
  p.temp    = NULL;                                                     \
  p.m_row = m;                                                          \
  p.m_col = n;                                                          \
  p.sizeof_T = sizeof(T);                                               \
  p.m_pdata = malloc(p.sizeof_T*m*n);                                   \
/*  memset(p.m_pdata,value,p.sizeof_T*m*n);   */                        \
  Matrix_init(p, value);                                                \  
} while(0)
                                                

#define Matrix_init_w_array(p, m, n, q) do {                            \
  Matrix_redim(p, m, n);                                                \
  long a, b;                                                            \
  for(a = 1; a <= p.m_row; a++){                                        \
    for(b = 1; b <= p.m_col; b++){                                      \
      Mat_v(p, a, b) =  q[(a-1)*n + (b-1)];                             \
    }                                                                   \
  }                                                                     \
} while(0)

/* A = delta_ij */
#define Matrix_eye(A, m) do {                                           \
  long I;                                                               \
  Matrix_redim(A, m, m);                                                \
  Matrix_init(A, 0.0);                                                  \
  for(I = 1; I <= m; I++)                                               \
    Mat_v(A, I, I) = 1.0;                                               \
                                                                        \
} while(0)

/* trA = tr_(A), trA = A_ii */
#define Matrix_trace(A,trA) do {                                        \
  long I;                                                               \
  trA = 0.0;                                                            \
  if(A.m_row != A.m_col || A.m_row ==0 || A.m_col ==0)                  \
    break;                                                              \ 
  for(I = 1; I <= A.m_row; I++)                                         \  
    trA += Mat_v(A, I, I);                                              \
} while(0)

#define Matrix_det(A, ddet) do {                                        \
                                                                        \
  if(A.m_row != A.m_col || A.m_row ==0 || A.m_col ==0)                  \
    break;                                                              \
                                                                        \
  if(A.m_row==1){                                                       \
    ddet = Mat_v(A, 1, 1);                                              \
    break;                                                              \
  }                                                                     \
                                                                        \
  if(A.m_row==2){                                                       \
    ddet = Mat_v(A, 1, 1)*Mat_v(A, 2, 2);                               \
    ddet -= Mat_v(A, 1, 2)*Mat_v(A, 2, 1);                              \
    break;                                                              \
  };                                                                    \
                                                                        \
  if(A.m_row==3){                                                       \
    ddet  = Mat_v(A, 1, 1)*Mat_v(A, 2, 2)*Mat_v(A, 3, 3);               \
    ddet += Mat_v(A, 1, 2)*Mat_v(A, 2, 3)*Mat_v(A, 3, 1);               \
    ddet += Mat_v(A, 1, 3)*Mat_v(A, 3, 2)*Mat_v(A, 2, 1);               \
    ddet -= Mat_v(A, 1, 3)*Mat_v(A, 2, 2)*Mat_v(A, 3, 1);               \
    ddet -= Mat_v(A, 1, 2)*Mat_v(A, 2, 1)*Mat_v(A, 3, 3);               \
    ddet -= Mat_v(A, 1, 1)*Mat_v(A, 3, 2)*Mat_v(A, 2, 3);               \
    break;                                                              \
  };                                                                    \
  if(A.m_row > 3){                                                      \
    printf("Matrix greater than [3x3] is not currently supported\n");   \
    break;                                                              \
  }                                                                     \
} while(0)                                                              

#define Matrix_inv(A, invA) do {                                        \
  double detA = 0.0;                                                    \
                                                                        \
  if(A.m_row != A.m_col || A.m_row ==0 || A.m_col ==0)                  \
    break;                                                              \
                                                                        \
  Matrix_redim(invA, A.m_row, A.m_col);                                 \
  if(A.m_row==1){                                                       \
    Mat_v(invA, 1, 1) = 1.0/Mat_v(A, 1, 1);                             \
    break;                                                              \
  }                                                                     \
                                                                        \
  Matrix_det(A, detA);                                                  \
  if(fabs(detA)<1.0e-15){                                               \
    printf("det(Matrix) = %f is close to zero\n", fabs(detA));          \
    break;                                                              \
  }                                                                     \
                                                                        \
  if(A.m_row==2){                                                       \
    Mat_v(invA, 1, 1) = Mat_v(A, 2, 2)/detA;                            \
    Mat_v(invA, 2, 2) = Mat_v(A, 1, 1)/detA;                            \
    Mat_v(invA, 1, 2) = -Mat_v(A, 1, 2)/detA;                           \
    Mat_v(invA, 2, 1) = -Mat_v(A, 2, 1)/detA;                           \
    break;                                                              \
  };                                                                    \
                                                                        \
  if(A.m_row==3){                                                       \
                                                                        \
    Mat_v(invA, 1, 1) = (Mat_v(A, 2, 2)*Mat_v(A, 3, 3)                  \
                       - Mat_v(A, 2, 3)*Mat_v(A, 3, 2))/detA;           \
    Mat_v(invA, 1, 2) = (Mat_v(A, 1, 3)*Mat_v(A, 3, 2)                  \
                       - Mat_v(A, 1, 2)*Mat_v(A, 3, 3))/detA;           \
    Mat_v(invA, 1, 3) = (Mat_v(A, 1, 2)*Mat_v(A, 2, 3)                  \
                       - Mat_v(A, 1, 3)*Mat_v(A, 2, 2))/detA;           \
    Mat_v(invA, 2, 1) = (Mat_v(A, 2, 3)*Mat_v(A, 3, 1)                  \
                       - Mat_v(A, 2, 1)*Mat_v(A, 3, 3))/detA;           \
    Mat_v(invA, 2, 2) = (Mat_v(A, 1, 1)*Mat_v(A, 3, 3)                  \
                       - Mat_v(A, 1, 3)*Mat_v(A, 3, 1))/detA;           \
    Mat_v(invA, 2, 3) = (Mat_v(A, 1, 3)*Mat_v(A, 2, 1)                  \
                       - Mat_v(A, 1, 1)*Mat_v(A, 2, 3))/detA;           \
    Mat_v(invA, 3, 1) = (Mat_v(A, 2, 1)*Mat_v(A, 3, 2)                  \
                       - Mat_v(A, 2, 2)*Mat_v(A, 3, 1))/detA;           \
    Mat_v(invA, 3, 2) = (Mat_v(A, 1, 2)*Mat_v(A, 3, 1)                  \
                       - Mat_v(A, 1, 1)*Mat_v(A, 3, 2))/detA;           \
    Mat_v(invA, 3, 3) = (Mat_v(A, 1, 1)*Mat_v(A, 2, 2)                  \
                       - Mat_v(A, 1, 2)*Mat_v(A, 2, 1))/detA;           \
    break;                                                              \
  };                                                                    \
  if(A.m_row > 3){                                                      \
    printf("Matrix greater than [3x3] is not currently supported\n");   \
    break;                                                              \
  }                                                                     \
} while(0)                                                              

#define Matrix_cleanup(A) do {                                          \
  A.m_row = 0;                                                          \
  A.m_col = 0;                                                          \
  if(A.m_pdata)                                                         \
    free(A.m_pdata);                                                    \
  A.m_pdata = NULL;                                                     \
  if(A.temp)                                                            \
    free(A.temp);                                                       \
  A.temp = NULL;                                                        \
} while(0)

#define Matrix_print(A) do {                                            \
  long a, b;                                                            \
  printf("[%ldx%ld] = \n", A.m_row, A.m_col);                           \
  for(a = 1; a <= A.m_row; a++){                                        \
    for(b = 1; b <= A.m_col; b++){                                      \
      printf("%e ", (double)Mat_v(A, a, b));                            \
    }                                                                   \
    printf("\n");                                                       \
  }                                                                     \
} while(0)

/* A = bB */
#define Matrix_AeqB(A, b, B) do {                                       \
  long I;                                                               \
  long m_row = B.m_row;                                                 \
  long m_col = B.m_col;                                                 \
  Matrix_redim(A, m_row, m_col);                                        \
                                                                        \
  for(I = 0; I < m_row*m_col; I++)                                      \
    A.m_pdata[I] = B.m_pdata[I]*b;                                      \
} while(0)
  
/* A = transpose(A) */
#define Matrix_trans(A) do {                                            \
  long a, b;                                                            \
  long m_row = A.m_row;                                                 \
  long m_col = A.m_col;                                                 \
  A.temp =  malloc(A.sizeof_T*m_row*m_col);                             \
                                                                        \
  for(a = 1; a <= m_row; a++){                                          \
    for(b = 1; b <= m_col; b++){                                        \
      A.temp[(b-1)*A.m_col+(a-1)] = Mat_v(A, a, b);                     \
    }                                                                   \
  }                                                                     \
  A.m_row = m_col;                                                      \
  A.m_col = m_row;                                                      \
  if(A.m_pdata)                                                         \
    free(A.m_pdata);                                                    \
  A.m_pdata = A.temp;                                                   \
  A.temp = NULL;                                                        \
} while(0)

/* C = aA + bB */
#define Matrix_AplusB(C, a, A, b, B) do {                               \
  long I, J;                                                            \
  long m_row = A.m_row;                                                 \
  long m_col = A.m_col;                                                 \
  if(A.m_row != B.m_row)                                                \
    break;                                                              \
  if(A.m_col != B.m_col)                                                \
    break;                                                              \
  Matrix_redim(C, m_row, m_col);                                        \
                                                                        \
  for(I = 0; I < m_row*m_col; I++)                                      \
    C.m_pdata[I] = A.m_pdata[I]*a + B.m_pdata[I]*b;                     \
} while(0)


/*C = aAxB + bC*/
// C[m,n] = a*A[m,k] x B[k,n] + b*C[m,n]
// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,a,A,k,B,n,b,C,n);
#define Matrix_AxB(C, a, b, A, AT, B, BT) do {                                                                                             \
  if(AT==0 && BT==0){                                                                                                                      \
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,A.m_row,B.m_col,A.m_col,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);\
  }                                                                                                                                        \
  if(AT==1 && BT==0){                                                                                                                      \
    cblas_dgemm(CblasRowMajor,CblasTrans,  CblasNoTrans,C.m_row,C.m_col,B.m_row,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);\
  }                                                                                                                                        \   
  if(AT==0 && BT==1){                                                                                                                      \
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,  C.m_row,C.m_col,A.m_col,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);\
  }                                                                                                                                        \
  if(AT==1 && BT==1){                                                                                                                      \
    cblas_dgemm(CblasRowMajor,CblasTrans,  CblasTrans,  C.m_row,C.m_col,A.m_row,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);\
  }                                                                                                                                        \
} while(0)

#endif
