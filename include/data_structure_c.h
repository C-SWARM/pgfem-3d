#ifndef _datastructure_c_H_
#define _datastructure_c_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mkl_cblas.h"
#include "utils.h"

#define Define_Matrix(T)                                                \
typedef struct Matrix_##T                                               \
{                                                                       \
  long m_row, m_col;                                                    \
  long sizeof_T;                                                        \
  T *m_pdata;                                                           \
  T *temp;                                                              \
} Matrix_##T                                                            \


#define Matrix(T) Matrix_##T
#define Mat_v(p, m, n) (p).m_pdata[((m)-1)*(p).m_col+((n)-1)]
#define Vec_v(p, m) (p).m_pdata[(m)-1]
#define Tns4_v(p, I,J,K,L) Vec_v(p,(I-1)*3*3*3+(J-1)*3*3+(K-1)*3+L)

#define Matrix_Tns4_mat_9x9(p) do{                                      \
  (p).m_row = 9;                                                        \
  (p).m_col = 9;                                                        \
} while(0)

#define Matrix_Mat2Vec(p) do{                                           \
  (p).m_row = (p).m_row*(p).m_col;                                      \
  (p).m_col = 1;                                                        \
} while(0)

#define Matrix_Vec2Mat(p, m, n) do{                                     \
  (p).m_row = (m);                                                      \
  (p).m_col = (n);                                                      \
} while(0)

#define Matrix_Tns4_eye(p) do{                                          \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Tns4_v(p, I,J,K,L) = 1.0*(I==K)*(J==L);                       \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

#define Matrix_Tns4_eye_bar(p) do{                                      \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Tns4_v(p, I,J,K,L) = 1.0*(I==L)*(J==K);                       \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

#define Matrix_construct(T, p) do {                                     \
  (p).m_pdata = NULL;                                                   \
  (p).temp    = NULL;                                                   \
  (p).m_row = 0;                                                        \
  (p).m_col = 0;                                                        \
  (p).sizeof_T = sizeof(T);                                             \
} while(0) 

#define Matrix_construct_redim(T, p, m, n) do {                         \
  (p).temp    = NULL;                                                   \
  (p).m_row = m;                                                        \
  (p).m_col = n;                                                        \
  (p).sizeof_T = sizeof(T);                                             \
  (p).m_pdata = malloc((p).sizeof_T*(m)*(n));                           \
} while(0)

#define Matrix_redim(p, m, n) do {                                      \
  Matrix_cleanup(p);                                                    \
  (p).m_pdata = malloc((p).sizeof_T*(m)*(n));                           \
  (p).temp    = NULL;                                                   \
  (p).m_row = m;                                                        \
  (p).m_col = n;                                                        \
} while(0)                                     

#define Matrix_init(p, value) do {                                      \
  long __a;                                                             \
  for(__a = 0; __a < (p).m_row*(p).m_col; __a++){                       \
      (p).m_pdata[__a] =  value;                                        \
  }                                                                     \
} while(0) 

#define Matrix_construct_init(T, p, m, n, value) do {                   \
  (p).temp    = NULL;                                                   \
  (p).m_row = m;                                                        \
  (p).m_col = n;                                                        \
  (p).sizeof_T = sizeof(T);                                             \
  (p).m_pdata = malloc((p).sizeof_T*(m)*(n));                           \
/*  memset(p.m_pdata,value,p.sizeof_T*m*n);   */                        \
  Matrix_init(p, value);                                                \
} while(0)
                                                

#define Matrix_init_w_array(p, m, n, q) do {                            \
  Matrix_redim(p, m, n);                                                \
  long __a, __b;                                                        \
  for(__a = 1; __a <= p.m_row; __a++){                                  \
    for(__b = 1; __b <= p.m_col; __b++){                                \
      Mat_v(p, __a, __b) =  q[(__a-1)*(n) + (__b-1)];                   \
    }                                                                   \
  }                                                                     \
} while(0)

/* A = delta_ij */
#define Matrix_eye(A, m) do {                                           \
  long __a;                                                             \
  Matrix_redim(A, m, m);                                                \
  Matrix_init(A, 0.0);                                                  \
  for(__a = 1; __a <= m; __a++)                                         \
    Mat_v(A, __a, __a) = 1.0;                                           \
                                                                        \
} while(0)

/* trA = tr_(A), trA = A_ii */
#define Matrix_trace(A,trA) do {                                        \
  long __a;                                                             \
  trA = 0.0;                                                            \
  if(A.m_row != A.m_col || A.m_row ==0 || A.m_col ==0)                  \
    break;                                                              \ 
  for(__a = 1; __a <= A.m_row; __a++)                                   \  
    trA += Mat_v(A, __a, __a);                                          \
} while(0)

#define Matrix_det(A, ddet) do {                                        \
                                                                        \
  if((A).m_row != (A).m_col || (A).m_row ==0 || (A).m_col ==0)          \
    break;                                                              \
                                                                        \
  if((A).m_row==1){                                                     \
    ddet = Mat_v(A, 1, 1);                                              \
    break;                                                              \
  }                                                                     \
                                                                        \
  if((A).m_row==2){                                                     \
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
  if((A).m_row > 3){                                                    \
    printf("Matrix greater than [3x3] is not currently supported\n");   \
    break;                                                              \
  }                                                                     \
} while(0)                                                              

#define Matrix_inv(A, invA) do {                                        \
  inverse((A).m_pdata,(A).m_row,(invA).m_pdata);                        \
} while(0)

#define Matrix_inv_no_use(A, invA) do {                                 \
  double detA = 0.0;                                                    \
                                                                        \
  if((A).m_row != (A).m_col || (A).m_row ==0 || (A).m_col ==0)          \
    break;                                                              \
                                                                        \
  Matrix_redim(invA, (A).m_row, (A).m_col);                             \
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
  if((A).m_row==2){                                                     \
    Mat_v(invA, 1, 1) =  Mat_v(A, 2, 2)/detA;                           \
    Mat_v(invA, 2, 2) =  Mat_v(A, 1, 1)/detA;                           \
    Mat_v(invA, 1, 2) = -Mat_v(A, 1, 2)/detA;                           \
    Mat_v(invA, 2, 1) = -Mat_v(A, 2, 1)/detA;                           \
    break;                                                              \
  };                                                                    \
                                                                        \
  if((A).m_row==3){                                                     \
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
  if((A).m_row > 3){                                                    \
    printf("Matrix greater than [3x3] is not currently supported\n");   \
    break;                                                              \
  }                                                                     \
} while(0)                                                              

#define Matrix_cleanup(A) do {                                          \
  (A).m_row = 0;                                                        \
  (A).m_col = 0;                                                        \
  if((A).m_pdata)                                                       \
    free((A).m_pdata);                                                  \
  (A).m_pdata = NULL;                                                   \
  if((A).temp)                                                          \
    free((A).temp);                                                     \
  (A).temp = NULL;                                                      \
} while(0)

#define Matrix_print(A) do {                                            \
  long __a, __b;                                                        \
  printf("[%ldx%ld] = \n", (A).m_row, (A).m_col);                       \
  for(__a = 1; __a <= (A).m_row; __a++){                                \
    for(__b = 1; __b <= (A).m_col; __b++){                              \
      printf("%e ", (double)Mat_v(A, __a, __b));                        \
    }                                                                   \
    printf("\n");                                                       \
  }                                                                     \
} while(0)

/* A = bB */
#define Matrix_AeqB(A, b, B) do {                                       \
  long __a;                                                             \
  long m_row = (B).m_row;                                               \
  long m_col = (B).m_col;                                               \
  Matrix_redim(A, m_row, m_col);                                        \
                                                                        \
  for(__a = 0; __a < m_row*m_col; __a++)                                \
    A.m_pdata[__a] = (B).m_pdata[__a]*(b);                              \
} while(0)
  
/* A = transpose(A) */
#define Matrix_trans(A) do {                                            \
  long __a, __b;                                                        \
  long m_row = (A).m_row;                                               \
  long m_col = (A).m_col;                                               \
  (A).temp =  malloc((A).sizeof_T*m_row*m_col);                         \
                                                                        \
  for(__a = 1; __a <= m_row; __a++){                                    \
    for(__b = 1; __b <= m_col; __b++){                                  \
      (A).temp[(__b-1)*(A).m_col+(__a-1)] = Mat_v(A, __a, __b);         \
    }                                                                   \
  }                                                                     \
  (A).m_row = m_col;                                                    \
  (A).m_col = m_row;                                                    \
  if((A).m_pdata)                                                       \
    free((A).m_pdata);                                                  \
  (A).m_pdata = (A).temp;                                               \
  (A).temp = NULL;                                                      \
} while(0)

/* C = aA + bB */
#define Matrix_AplusB(C, a, A, b, B) do {                               \
  long __I;                                                             \
  long m_row = (A).m_row;                                               \
  long m_col = (A).m_col;                                               \
  if((A).m_row != (B).m_row)                                            \
    break;                                                              \
  if((A).m_col != (B).m_col)                                            \
    break;                                                              \
                                                                        \
  for(__I = 0; __I < m_row*m_col; __I++)                                \
    (C).m_pdata[__I] = (A).m_pdata[__I]*(a) + (B).m_pdata[__I]*(b);     \
} while(0)

#define Matrix_AOxB(C, A, B) do {                                       \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Tns4_v(C, I,J,K,L) = Mat_v(A,I,J)*Mat_v(B,K,L);               \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

#define Matrix_Tns4_dd_Tns2(C, A, B) do {                               \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      Mat_v(C,I,J)=0.0;                                                 \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Mat_v(C,I,J) += Tns4_v(A,I,J,K,L)*Mat_v(B,K,L);               \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

/*C = aAxB + bC*/
// C[m,n] = a*A[m,k] x B[k,n] + b*C[m,n]
// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,a,A,k,B,n,b,C,n);
#define Matrix_AxB(C, a, b, A, AT, B, BT) do {                                                                                                               \
  if(AT==0 && BT==0){                                                                                                                                        \
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(A).m_row,(B).m_col,(A).m_col,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
  if(AT==1 && BT==0){                                                                                                                                        \
    cblas_dgemm(CblasRowMajor,CblasTrans,  CblasNoTrans,(C).m_row,(C).m_col,(B).m_row,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \   
  if(AT==0 && BT==1){                                                                                                                                        \
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,  (C).m_row,(C).m_col,(A).m_col,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
  if(AT==1 && BT==1){                                                                                                                                        \
    cblas_dgemm(CblasRowMajor,CblasTrans,  CblasTrans,  (C).m_row,(C).m_col,(A).m_row,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
} while(0)

#endif
