#ifndef _datastructure_c_H_
#define _datastructure_c_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "index_macros.h"
#include "mkl_cblas.h"
#include "utils.h"
#include <algorithm>
#include <numeric>

// Tolerate the legacy macros that are too much trouble to replace.
#define Matrix_construct(T, p) (p).construct()
#define Matrix_construct_redim(T, p, m, n) (p).construct(m, n)
#define Matrix_construct_init(T, p, m, n, value) (p).construct(m, n, value)

template <class T>
class Matrix {
 public:
  // static constexpr size_t sizeof_T = sizeof(T);
  size_t m_row = 0;
  size_t m_col = 0;
  T *m_pdata = nullptr;

  Matrix& operator=(const Matrix& rhs) {
    assert(m_row == rhs.m_row and m_col == rhs.m_col);
    std::copy_n(rhs.m_pdata, m_row * m_col, m_pdata);
    return *this;
  }

  template <class U>
  Matrix<T>& operator=(const Matrix<U>& rhs) {
    assert(m_row == rhs.m_row and m_col == rhs.m_col);
    std::copy_n(rhs.m_pdata, m_row * m_col, m_pdata);
    return *this;
  }

  void construct() {
    assert(!m_pdata);
    m_row = 0;
    m_col = 0;
    m_pdata = nullptr;
  }

  void construct(size_t m, size_t n) {
    if (m_pdata) {
      delete [] m_pdata;
    }
    m_row = m;
    m_col = n;
    m_pdata = new T[m * n];
  }

  template <class U>
  void construct(size_t m, size_t n, U&& val) {
    construct(m, n);
    std::fill_n(m_pdata, m * n, std::forward<U>(val));
  }

  template <class U>
  void init(const U* array) {
    assert(m_pdata);
    std::copy_n(array, m_row * m_col, m_pdata);
  }

  template <class U>
  void init(U&& val) {
    assert(m_pdata);
    std::fill_n(m_pdata, m_row * m_col, val);
  }

  void init() {
    init(T{});
  }

  void cleanup() {
    if (m_pdata) {
      delete [] m_pdata;
    }
    m_row = 0;
    m_col = 0;
    m_pdata = nullptr;
  }

  void redim(size_t m, size_t n) {
    cleanup();
    construct(m, n);
  }

  void check_null_and_redim(size_t m, size_t n) {
    if ((m_pdata == nullptr) || (m_row != m) || (m_col != n)) {
      redim(m, n);
    }
  }

  void trans() {
    T* temp = new T[m_col * m_row];
    for (unsigned i = 0; i < m_row; ++i) {
      for (unsigned j = 0; j < m_col; ++j) {
        // temp[j][i] = m_pdata[i][j];
        temp[j*m_row + i] = m_pdata[i*m_col +j];
      }
    }
    std::swap(m_col, m_row);
    std::swap(m_pdata, temp);
    delete [] temp;
  }

  T& Vec_v(size_t m) {
    return m_pdata[m-1];
  }

  const T& Vec_v(size_t m) const {
    return m_pdata[m-1];
  }

  T& Mat_v(size_t m, size_t n) {
    return m_pdata[idx_2_gen(m - 1, n - 1, m_row, m_col)];
  }

  const T& Mat_v(size_t m, size_t n) const {
    return m_pdata[idx_2_gen(m - 1, n - 1, m_row, m_col)];
  }

  T& Tns4_v(size_t I, size_t J, size_t K, size_t L) {
    return m_pdata[idx_4((I) - 1, (J) - 1, (K) - 1, (L) - 1)];
  }

  const T& Tns4_v(size_t I, size_t J, size_t K, size_t L) const {
    return m_pdata[idx_4((I) - 1, (J) - 1, (K) - 1, (L) - 1)];
  }

  void Tns4_mat_9x9() {
    m_row = m_col = 9;
  }

  void Mat2Vec() {
    m_row = m_row * m_col;
    m_col = 1;
  }

  void Vec2Mat(size_t m, size_t n) {
    m_row = m;
    m_col = n;
  }

  void Tns4_eye() {
    assert(m_row == 9);
    assert(m_col == 9);
    init();
    for(int I=1; I<=3; I++)
      for(int J=1; J<=3; J++)
        for(int K=1; K<=3; K++)
          for(int L=1; L<=3; L++)
            Tns4_v(I,J,K,L) = 1.0*(I==K)*(J==L);
  }

  void Tns4_II() {
    assert(m_row == 9);
    assert(m_col == 9);
    init();
    for(int I=1; I<=3; I++)
      for(int J=1; J<=3; J++)
        for(int K=1; K<=3; K++)
          for(int L=1; L<=3; L++)
            Tns4_v(I,J,K,L) = 0.5*((I==K)*(J==L)+(I==L)*(J==K));
  }

  void Tns4_eye_bar() {
    assert(m_row == 9);
    assert(m_col == 9);
    init();
    for(int I=1; I<=3; I++)
      for(int J=1; J<=3; J++)
        for(int K=1; K<=3; K++)
          for(int L=1; L<=3; L++)
            Tns4_v(I,J,K,L) = 1.0*(I==L)*(J==K);
  }

  void eye() {
    assert(m_row == m_col);
    init();
    for (unsigned i = 1; i <= m_row; ++i) {
      Mat_v(i,i) = T{1};
    }
  }

  T trace() const {
    T out{};
    if(m_row != m_col || m_row ==0 || m_col ==0) {
      return out;
    }
    for (unsigned i = 1; i <= m_row; ++i) {
      out += Mat_v(i,i);
    }
    return out;
  }

  T det() const {
    if (m_row != m_col || m_row ==0 || m_col ==0)
      return T{0};

    if(m_row==1){
      return Mat_v(1, 1);
    }

    if(m_row==2){
      return Mat_v(1, 1) * Mat_v(2, 2) -
             Mat_v(1, 2) * Mat_v(2, 1);
    }

    if(m_row==3){
      return Mat_v(1, 1) * Mat_v(2, 2) * Mat_v(3, 3) +
             Mat_v(1, 2) * Mat_v(2, 3) * Mat_v(3, 1) +
             Mat_v(1, 3) * Mat_v(3, 2) * Mat_v(2, 1) -
             Mat_v(1, 3) * Mat_v(2, 2) * Mat_v(3, 1) -
             Mat_v(1, 2) * Mat_v(2, 1) * Mat_v(3, 3) -
             Mat_v(1, 1) * Mat_v(3, 2) * Mat_v(2, 3);
    }

    printf("Matrix greater than [3x3] is not currently supported\n");
    return T{0};
  }

  int inverse(Matrix<T>& out) const {
    return ::inverse(m_pdata, m_row, out.m_pdata);
  }

  void inverse_no_use(Matrix<T>& out) const {
    if(m_row != m_col || m_row ==0 || m_col ==0)
      return;

    out.check_null_and_redim(m_row, m_col);
    if(m_row==1){
      out.Mat_v(1, 1) = 1.0/Mat_v(1, 1);
      return;
    }

    auto detA = det();
    if(fabs(detA)<1.0e-15){
      printf("det(Matrix) = %f is close to zero\n", fabs(detA));
      return;
    }

    auto invD = 1/detA;

    if(m_row==2){
      out.Mat_v(1, 1) =  Mat_v(2, 2) * invD;
      out.Mat_v(2, 2) =  Mat_v(1, 1) * invD;
      out.Mat_v(1, 2) = -Mat_v(1, 2) * invD;
      out.Mat_v(2, 1) = -Mat_v(2, 1) * invD;
      return;
    }

    if(m_row==3){
      out.Mat_v(1, 1) = (Mat_v(2, 2)*Mat_v(3, 3) - Mat_v(2, 3)*Mat_v(3, 2))*invD;
      out.Mat_v(1, 2) = (Mat_v(1, 3)*Mat_v(3, 2) - Mat_v(1, 2)*Mat_v(3, 3))*invD;
      out.Mat_v(1, 3) = (Mat_v(1, 2)*Mat_v(2, 3) - Mat_v(1, 3)*Mat_v(2, 2))*invD;
      out.Mat_v(2, 1) = (Mat_v(2, 3)*Mat_v(3, 1) - Mat_v(2, 1)*Mat_v(3, 3))*invD;
      out.Mat_v(2, 2) = (Mat_v(1, 1)*Mat_v(3, 3) - Mat_v(1, 3)*Mat_v(3, 1))*invD;
      out.Mat_v(2, 3) = (Mat_v(1, 3)*Mat_v(2, 1) - Mat_v(1, 1)*Mat_v(2, 3))*invD;
      out.Mat_v(3, 1) = (Mat_v(2, 1)*Mat_v(3, 2) - Mat_v(2, 2)*Mat_v(3, 1))*invD;
      out.Mat_v(3, 2) = (Mat_v(1, 2)*Mat_v(3, 1) - Mat_v(1, 1)*Mat_v(3, 2))*invD;
      out.Mat_v(3, 3) = (Mat_v(1, 1)*Mat_v(2, 2) - Mat_v(1, 2)*Mat_v(2, 1))*invD;
      return;
    }

    printf("Matrix greater than [3x3] is not currently supported\n");
    return;
  }

  void print() const {
    printf("[%ldx%ld] = \n", m_row, m_col);
    for(unsigned i = 1; i <= m_row; ++i){
      for(unsigned j = 1; j <= m_col; ++j){
        printf("%e ", Mat_v(i, j));
      }
      printf("\n");
    }
  }

  void print_name(const char name[]) const {
    printf("%s = [\n", name);
    for(unsigned a__A = 1; a__A <= m_row; a__A++){
      for(unsigned b__B = 1; b__B <= m_col; b__B++){
        printf("%e ", Mat_v(a__A, b__B));
      }
      if(a__A==m_row)
        printf("];\n");
      else
        printf("\n");
    }
  }

  /* A = bB */
  template <class U, class V>
  void eqB(U&& b, const Matrix<V>& rhs) {
    check_null_and_redim(rhs.m_row, rhs.m_col);
    for (unsigned i = 0, e = m_row * m_col; i < e; ++i) {
      m_pdata[i] = b * rhs.m_pdata[i];
    }
  }

  /* A = bBT */
  template <class U, class V>
  void eqBT(U&& b, const Matrix<V>& rhs) {
    check_null_and_redim(rhs.m_row, rhs.m_col);
    for(unsigned MaTtEmPVar_a = 1; MaTtEmPVar_a <= m_row; MaTtEmPVar_a++)
      for(unsigned MaTtEmPVar_b = 1; MaTtEmPVar_b <= m_col; MaTtEmPVar_b++)
        Mat_v(MaTtEmPVar_b, MaTtEmPVar_a) = rhs.Mat_v(MaTtEmPVar_a, MaTtEmPVar_b)*b;
  }

  /* C = aA + bB */
  template <class N, class M, class O, class P>
  void eqAPlusB(N&& a, const Matrix<M>& A, O&& b, const Matrix<P>& B) {
    if(A.m_row != B.m_row || A.m_col != B.m_col)
      return;
    check_null_and_redim(A.m_row, A.m_col);
    for(unsigned i = 0, e = m_row * m_col; i < e; ++i)
      m_pdata[i] = A.m_pdata[i]*a + B.m_pdata[i]*b;
  }

  template <class U, class V>
  void eqAOxB(const Matrix<U>& A, const Matrix<V>& B) {
    assert(m_col == 9);
    assert(m_row == 9);
    for(int I=1; I<=3; I++)
      for(int J=1; J<=3; J++)
        for(int K=1; K<=3; K++)
          for(int L=1; L<=3; L++)
            Tns4_v(I,J,K,L) = A.Mat_v(I,J)*B.Mat_v(K,L);
  }

  template <class U, class V>
  void eqTns4_dd_Tns2(const Matrix<U>& A, const Matrix<V>& B) {
    init(0.0);
    for(int K=1; K<=3; K++)
      for(int L=1; L<=3; L++)
        if(fabs(B.Mat_v(K,L))>=1.0e-15)
          for(int I=1; I<=3; I++)
            for(int J=1; J<=3; J++)
              Mat_v(I,J) += A.Tns4_v(I,J,K,L)*B.Mat_v(K,L);
  }

  template <class U, class V>
  void eqTns2_dd_Tns4(const Matrix<U>& A, const Matrix<V>& B) {
    init(0.0);
    for(int I=1; I<=3; I++)
      for(int J=1; J<=3; J++)
        if(A.Mat_v(I,J)>=1.0e-15)
          for(int K=1; K<=3; K++)
            for(int L=1; L<=3; L++)
              Mat_v(K,L) += A.Mat_v(I,J)*B.Tns4_v(I,J,K,L);
  }

  template <class U>
  T ddot(const Matrix<U>& B) const {
    return cblas_ddot(m_row * m_col, m_pdata, 1, B.m_pdata, 1);
  }

  template <class U>
  T ddot_no_use(const Matrix<U>& B) const {
    return std::inner_product(m_pdata, m_pdata + m_row * m_col, B.m_pdata, 0.0);
  }

  template <class U, class V, class W, class X>
  void eqGEMM(U&& a, V&& b, const Matrix<W>& A, bool AT,
                            const Matrix<X>& B, bool BT) {
    if(AT==0 && BT==0){
      check_null_and_redim(A.m_row, B.m_col);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.m_row, B.m_col,
          A.m_col, a, A.m_pdata, A.m_col, B.m_pdata, B.m_col, b, m_pdata, m_col);
    }
    if(AT==1 && BT==0){
      check_null_and_redim(A.m_col, B.m_col);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m_row, m_col,
          B.m_row, a, A.m_pdata, A.m_col, B.m_pdata, B.m_col, b, m_pdata,
          m_col);
    }
    if(AT==0 && BT==1){
      check_null_and_redim(A.m_row, B.m_row);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m_row, m_col,
          A.m_col, a, A.m_pdata, A.m_col, B.m_pdata, B.m_col, b, m_pdata,
          m_col);
    }
    if(AT==1 && BT==1){
      check_null_and_redim(A.m_col, B.m_row);
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m_row, m_col, A.m_row,
          a, A.m_pdata, A.m_col, B.m_pdata, B.m_col, b, m_pdata, m_col);
    }
  }

  // D = a*A*B*C + b*D, a and b are scalar
  template <class U, class V, class W, class X, class Y>
  void eqTns2_AxBxC(U&& a, V&& b, const Matrix<W>& A, const Matrix<X>& B,
                             const Matrix<Y>& C)  {
    check_null_and_redim(3,3);
    Mat_v(1,1) = Mat_v(1,1)*b + a*(C.Mat_v(1,1)*(A.Mat_v(1,1)*B.Mat_v(1,1) + A.Mat_v(1,2)*B.Mat_v(2,1) + A.Mat_v(1,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,1)*(A.Mat_v(1,1)*B.Mat_v(1,2) + A.Mat_v(1,2)*B.Mat_v(2,2) + A.Mat_v(1,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,1)*(A.Mat_v(1,1)*B.Mat_v(1,3) + A.Mat_v(1,2)*B.Mat_v(2,3) + A.Mat_v(1,3)*B.Mat_v(3,3)));
    Mat_v(1,2) = Mat_v(1,2)*b + a*(C.Mat_v(1,2)*(A.Mat_v(1,1)*B.Mat_v(1,1) + A.Mat_v(1,2)*B.Mat_v(2,1) + A.Mat_v(1,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,2)*(A.Mat_v(1,1)*B.Mat_v(1,2) + A.Mat_v(1,2)*B.Mat_v(2,2) + A.Mat_v(1,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,2)*(A.Mat_v(1,1)*B.Mat_v(1,3) + A.Mat_v(1,2)*B.Mat_v(2,3) + A.Mat_v(1,3)*B.Mat_v(3,3)));
    Mat_v(1,3) = Mat_v(1,3)*b + a*(C.Mat_v(1,3)*(A.Mat_v(1,1)*B.Mat_v(1,1) + A.Mat_v(1,2)*B.Mat_v(2,1) + A.Mat_v(1,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,3)*(A.Mat_v(1,1)*B.Mat_v(1,2) + A.Mat_v(1,2)*B.Mat_v(2,2) + A.Mat_v(1,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,3)*(A.Mat_v(1,1)*B.Mat_v(1,3) + A.Mat_v(1,2)*B.Mat_v(2,3) + A.Mat_v(1,3)*B.Mat_v(3,3)));
    Mat_v(2,1) = Mat_v(2,1)*b + a*(C.Mat_v(1,1)*(A.Mat_v(2,1)*B.Mat_v(1,1) + A.Mat_v(2,2)*B.Mat_v(2,1) + A.Mat_v(2,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,1)*(A.Mat_v(2,1)*B.Mat_v(1,2) + A.Mat_v(2,2)*B.Mat_v(2,2) + A.Mat_v(2,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,1)*(A.Mat_v(2,1)*B.Mat_v(1,3) + A.Mat_v(2,2)*B.Mat_v(2,3) + A.Mat_v(2,3)*B.Mat_v(3,3)));
    Mat_v(2,2) = Mat_v(2,2)*b + a*(C.Mat_v(1,2)*(A.Mat_v(2,1)*B.Mat_v(1,1) + A.Mat_v(2,2)*B.Mat_v(2,1) + A.Mat_v(2,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,2)*(A.Mat_v(2,1)*B.Mat_v(1,2) + A.Mat_v(2,2)*B.Mat_v(2,2) + A.Mat_v(2,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,2)*(A.Mat_v(2,1)*B.Mat_v(1,3) + A.Mat_v(2,2)*B.Mat_v(2,3) + A.Mat_v(2,3)*B.Mat_v(3,3)));
    Mat_v(2,3) = Mat_v(2,3)*b + a*(C.Mat_v(1,3)*(A.Mat_v(2,1)*B.Mat_v(1,1) + A.Mat_v(2,2)*B.Mat_v(2,1) + A.Mat_v(2,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,3)*(A.Mat_v(2,1)*B.Mat_v(1,2) + A.Mat_v(2,2)*B.Mat_v(2,2) + A.Mat_v(2,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,3)*(A.Mat_v(2,1)*B.Mat_v(1,3) + A.Mat_v(2,2)*B.Mat_v(2,3) + A.Mat_v(2,3)*B.Mat_v(3,3)));
    Mat_v(3,1) = Mat_v(3,1)*b + a*(C.Mat_v(1,1)*(A.Mat_v(3,1)*B.Mat_v(1,1) + A.Mat_v(3,2)*B.Mat_v(2,1) + A.Mat_v(3,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,1)*(A.Mat_v(3,1)*B.Mat_v(1,2) + A.Mat_v(3,2)*B.Mat_v(2,2) + A.Mat_v(3,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,1)*(A.Mat_v(3,1)*B.Mat_v(1,3) + A.Mat_v(3,2)*B.Mat_v(2,3) + A.Mat_v(3,3)*B.Mat_v(3,3)));
    Mat_v(3,2) = Mat_v(3,2)*b + a*(C.Mat_v(1,2)*(A.Mat_v(3,1)*B.Mat_v(1,1) + A.Mat_v(3,2)*B.Mat_v(2,1) + A.Mat_v(3,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,2)*(A.Mat_v(3,1)*B.Mat_v(1,2) + A.Mat_v(3,2)*B.Mat_v(2,2) + A.Mat_v(3,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,2)*(A.Mat_v(3,1)*B.Mat_v(1,3) + A.Mat_v(3,2)*B.Mat_v(2,3) + A.Mat_v(3,3)*B.Mat_v(3,3)));
    Mat_v(3,3) = Mat_v(3,3)*b + a*(C.Mat_v(1,3)*(A.Mat_v(3,1)*B.Mat_v(1,1) + A.Mat_v(3,2)*B.Mat_v(2,1) + A.Mat_v(3,3)*B.Mat_v(3,1))
                                   + C.Mat_v(2,3)*(A.Mat_v(3,1)*B.Mat_v(1,2) + A.Mat_v(3,2)*B.Mat_v(2,2) + A.Mat_v(3,3)*B.Mat_v(3,2))
                                   + C.Mat_v(3,3)*(A.Mat_v(3,1)*B.Mat_v(1,3) + A.Mat_v(3,2)*B.Mat_v(2,3) + A.Mat_v(3,3)*B.Mat_v(3,3)));
  }
};

template <class T>
using Vector = Matrix<T>;

/**
 * Indexing functions. !!NOTE!! Indexing starts from 1
 *
 * Provides direct read/write access to data.
 */
template <class T>
static inline T&
Vec_v(Matrix<T>& p, size_t m) {
  return p.Vec_v(m);
}

template <class T>
static inline const T&
Vec_v(const Matrix<T>& p, size_t m) {
  return p.Vec_v(m);
}

template <class T>
static inline T&
Mat_v(Matrix<T>& p, size_t m, size_t n) {
  return p.Mat_v(m, n);
}

template <class T>
static inline const T&
Mat_v(const Matrix<T>& p, size_t m, size_t n) {
  return p.Mat_v(m, n);
}

template <class T>
static inline const T&
Tns4_v(const Matrix<T>& p, size_t I, size_t J, size_t K, size_t L) {
  return p.Tns4_v(I, J, K, L);
}

template <class T>
static inline T&
Tns4_v(Matrix<T>& p, size_t I, size_t J, size_t K, size_t L) {
  return p.Tns4_v(I, J, K, L);
}

template <class T>
static inline void
Matrix_Tns4_mat_9x9(Matrix<T>& p) {
  p.Tns4_mat_9x9();
}

template <class T>
static inline void
Matrix_Mat2Vec(Matrix<T>& p) {
  p.Mat2Vec();
}

/**
 * Convert a (row) Vector to a [m x n] Matrix
 */
template <class T>
static inline void
Matrix_Vec2Mat(Matrix<T>& p, size_t m, size_t n) {
  p.Mat2Vec();
}

/**
 * Compute the 4th order identity tensor
 */
template <class T>
static inline void
Matrix_Tns4_eye(Matrix<T>& p) {
  p.Tns4_eye();
}

/**
 * Comptue the 4th order symmetric identity tensor
 */
template <class T>
static inline void
Matrix_Tns4_II(Matrix<T>& p) {
  p.Tns4_II();
}

/**
 * See Matrix_Tns4_eye, transpose k,l
 */
template <class T>
static inline void
Matrix_Tns4_eye_bar(Matrix<T>& p) {
  p.Tns4_eye_bar();
}

template <class T>
static inline void
Matrix_redim(Matrix<T>& p, size_t m, size_t n) {
  p.redim(m, n);
}

template <class T>
static inline void
Matrix_check_null_and_redim(Matrix<T>& p, size_t m, size_t n) {
  p.check_null_and_redim(m,n);
}

template <class T, class U>
static inline void
Matrix_init(Matrix<T>& p, U&& value) {
  p.init(std::forward<U>(value));
}

template <class T, class U>
static inline void
Matrix_init_w_array(Matrix<T>& p, size_t m, size_t n, const U* q) {
  p.check_null_and_redim(m, n);
  p.init(q);
}

/* A = delta_ij */
template <class T>
static inline void
Matrix_eye(Matrix<T>& A, size_t m) {
  Matrix_check_null_and_redim(A, m, m);
  A.eye();
}

// A_symm = symmetric(A)
// Matrix_symmetric(A,A_symm)
template <class T>
static inline void
Matrix_symmetric(Matrix<T>& A, Matrix<T>& A_symm) {
  symmetric_part(A_symm.m_pdata,A.m_pdata,A.m_row);
}

/* trA = tr_(A), trA = A_ii */
template <class T, class U>
static inline void
Matrix_trace(const Matrix<T>& A, U& trA) {
  trA = A.trace();
}

template <class T, class U>
static inline void
Matrix_det(const Matrix<T>& A, U& ddet) {
  ddet = A.det();
}

template <class T>
static inline void
Matrix_inv_check_err(const Matrix<T>& A, Matrix<T>& invA, int& info) {
  info = A.inverse(invA);
}

template <class T>
static inline void
Matrix_inv(const Matrix<T>& A, Matrix<T>& invA) {
  A.inverse(invA);
}

template <class T>
static inline void
Matrix_inv_no_use(const Matrix<T>& A, Matrix<T>& invA) {
  A.inverse_no_use(invA);
}

template <class T>
static inline void
Matrix_cleanup(Matrix<T>& A) {
  A.cleanup();
}

template <class T>
static inline void
Matrix_print(const Matrix<T>& A) {
  A.print();
}

template <class T>
static inline void
Matrix_print_name(const Matrix<T>& A, const char name[]) {
  A.print(name);
}

/* A <- B */
template <class T, class U>
static inline void
Matrix_copy(Matrix<T>& A, const Matrix<U>& B) {
  A = B;
}

/* A = bB */
template <class T, class U, class V>
static inline void
Matrix_AeqB(Matrix<T>& A, U&& b, const Matrix<V>& B) {
  A.eqB(std::forward<U>(b), B);
}

/* A = transpose(A) */
template <class T>
static inline void
Matrix_trans(Matrix<T>& A) {
  A.trans();
}

/* A = bBT */
template <class T, class U, class V>
static inline void
Matrix_AeqBT(Matrix<T>& A, U&& b, const Matrix<V>& B) {
  A.eqBT(std::forward<U>(b), B);
}

/* C = aA + bB */
template <class T, class U, class V, class W, class X>
static inline void
Matrix_AplusB(Matrix<T>& C, U&& a, const Matrix<V>& A,
                            W&& b, const Matrix<X>& B) {
  C.eqAPlusB(std::forward<U>(a), A, std::forward<W>(b), B);
}

template <class T, class U, class V>
static inline void
Matrix_AOxB(Matrix<T>& C, const Matrix<U>& A, const Matrix<V>& B) {
  C.eqAOxB(A, B);
}

// C = A:B
template <class T, class U, class V>
static inline void
Matrix_Tns4_dd_Tns2(Matrix<T>& C, const Matrix<U>& A, const Matrix<V>& B) {
  C.eqTns4_dd_Tns2(A, B);
}

// C = A:B
template <class T, class U, class V>
static inline void
Matrix_Tns2_dd_Tns4(Matrix<T>& C, const Matrix<U>& A, const Matrix<V>& B) {
  C.eqTns2_dd_Tns4(A, B);
}

template <class T, class U, class V>
static inline void
Matrix_ddot(const Matrix<T>& A, const Matrix<U>& B, V& out) {
  out = A.ddot(B);
}

template <class T, class U, class V>
static inline void
Matrix_ddot_no_use(const Matrix<T>& A, const Matrix<U>& B, V& out) {
  out = A.ddot_no_use(B);
}

/*C = aAxB + bC*/
// C[m,n] = a*A[m,k] x B[k,n] + b*C[m,n]
// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,a,A,k,B,n,b,C,n);
template <class T, class U, class V, class W, class X>
static inline void
Matrix_AxB(Matrix<T>& C, U&& a, V&& b, const Matrix<W>& A, bool AT,
                                       const Matrix<X>& B, bool BT) {
  C.eqGEMM(std::forward<U>(a), std::forward<V>(b), A, AT, B, BT);
}

// D = a*A*B*C + b*D, a and b are scalar
template <class T, class U, class V, class W, class X, class Y>
static inline void
Matrix_Tns2_AxBxC(Matrix<T>& D, U&& a, V&& b, const Matrix<W>& A,
                                              const Matrix<X>& B,
                                              const Matrix<Y>& C)  {
  D.eqTns2_AxBxC(std::forward<U>(a), std::forward<V>(b), A, B, C);
}

#endif
