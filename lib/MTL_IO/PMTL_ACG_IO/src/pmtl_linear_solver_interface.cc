#include<stdio.h>
#include "pmtl_linear_solver_interface.h"

#include <iostream>
#include <boost/mpi.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

using namespace mtl;
namespace mpi = boost::mpi;

#define mat matrix

#define Matrix(T)  mtl::dense2D<T>
#define Vector(T)  mtl::dense_vector<T> 
#define CMatrix(T) mtl::compressed2D<T>

#define Distribute(T) matrix::distributed<T>
#define DistributeV(T) vector::distributed<T>

typedef struct MTL_SOLVER_INTF {
  CMatrix(double) A;
  Distribute(CMatrix(double)) AA;
  Vector(double) x;
  Vector(double) b;  
  DistributeV(Vector(double)) bb;
  DistributeV(Vector(double)) xx;
  itl::pc::identity<CMatrix(double)> *P;
} MTL_SOLVER_INTF;

template <typename T>
int Matrix_print(T &A)
{
  int err = 0;
  printf("ans = [\n");
  for(int a=0; a<A.num_rows(); a++)
  {
    for(int b=0; b<A.num_cols(); b++)
      printf("%e ", A[a][b]);
    
    if(a==A.num_rows()-1)
      printf("];\n");
    else
      printf("\n");  
  }  
  return err;
}

void insert(Matrix(double) &m, Matrix(double) &input,int index1, int index2)
{
 if(input.num_cols()==0 || input.num_rows()==0)
   return;
   
 mat::inserter<Matrix(double)> ins(m);
 
 Vector(int) v1(input.num_rows());
 Vector(int) v2(input.num_cols());
 for (int i=0;i<input.num_rows();i++)
     v1[i] = index1 + i;
 
 for (int i=0;i<input.num_cols();i++) 
     v2[i] = index2 + i;
 
 ins << element_matrix(input,v1,v2);
}

int construct_solver(void **m)
{   
  int err = 0;
  MTL_SOLVER_INTF *temp = (MTL_SOLVER_INTF *) malloc(sizeof(MTL_SOLVER_INTF));
  if(temp == NULL)
  {
    printf("construct_solver: Out of memory!\n");
    exit(1);
  }
  
  *m = temp;

  return err;
}

int destruct_solver(void **m)
{   
  int err = 0;
  MTL_SOLVER_INTF *temp = (MTL_SOLVER_INTF *) *m;
  
  *m = NULL;
  
  free(temp);
  return err;
}

int initialize_linear_system(void *m, int N)
{ 
  int err = 0;  

  mpi::communicator comm;
  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;

  intf->A.change_dim(N,N);

  typedef mtl::matrix::distributed<mtl::compressed2D<double> >  matrix_type;
         matrix_type PA(N,N);

  typedef mtl::vector::distributed<mtl::dense_vector<double> >  vector_type;
  	  	 vector_type  Pb(N,0.0);

  PA = 0.0;
  intf->AA = PA;
  intf->b.change_dim(N);
  intf->x.change_dim(N); 
      
//   intf->AA = 0.0;
  intf->b = 0.0;
  intf->x = 0.0;
  intf->bb = Pb;
  intf->xx = Pb;
  return err;
}

int update_linear_system_A_IJ(void *m, int *I, int IN,
                                       int *J, int JN, double *values, MPI_Comm comm_1)
{   
  int err = 0;  

  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;
  boost::mpi::communicator comm(communicator(intf->AA));
  
  typedef Collection<CMatrix(double)>::value_type value_type;
  {
  mat::inserter<CMatrix(double), update_plus<value_type> > ins(intf->A, 3);

   
  Matrix(double) temp(IN,JN);
  Vector(double) vI(IN), vJ(JN);
  
  double *A = temp.address_data();

  for(int a=0; a<IN*JN; a++)
    A[a] = values[a];
    
  for(int a=0; a<IN; a++)
    vI[a] = I[a];
    
  for(int a=0; a<JN; a++)
    vJ[a] = J[a];   
  
  ins << element_matrix(temp, vI, vJ);
}
 // std :: cout << intf->b <<"\n";

  typedef mtl::matrix::distributed<mtl::compressed2D<double> >  matrix_type;
  typedef mtl::vector::distributed<mtl::dense_vector<double> >  vector_type;
  matrix_type PA(IN,JN);
  vector_type Pb(IN,0.0);
  {
  mtl::matrix::inserter<matrix_type, mtl::operations::update_plus<double> > insA(PA);
  mtl::vector::inserter<vector_type, mtl::operations::update_plus<double> > insV(Pb);

  for (int i=0;i<IN;i++){
	  insV[i] << intf->b[i];
     for (int j=0;j<JN;j++) {
       insA[i][j] << intf->A[i][j];
     }
   }
  }



  intf->AA=PA;
  intf->bb=Pb;
  //std :: cout << intf->AA <<"\n";
 // std :: cout << intf->bb <<"\n";
  return err;
}

int update_linear_system_A_IJ_(void *m, int *I, int IN,
                                       int *J, int JN, double *values)
{   
  int err = 0;  

  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;
  
  typedef Collection<CMatrix(double)>::value_type value_type;
  mat::inserter<CMatrix(double), update_plus<value_type> > ins(intf->A, 3);
    
  Matrix(double) temp(IN,JN);
  Vector(double) vI(IN), vJ(JN);
  
  double *A = temp.address_data();

  for(int a=0; a<IN*JN; a++)
    A[a] = values[a];
    
  for(int a=0; a<IN; a++)
    vI[a] = I[a];
    
  for(int a=0; a<JN; a++)
    vJ[a] = J[a];   
  
  ins << element_matrix(temp, vI, vJ);
  return err;
}

int set_linear_system_b(void *m, double *values, int N)
{   
  int err = 0;  

  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;

  double *b = (intf->b).address_data();
        
  for(int a=0; a<N; a++)
   b[a] = values[a];
  
  return err;
}

int set_solver_pc(void *m, int type)
{   
  int err = 0;

  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;  
  itl::pc::identity<CMatrix(double)> P(intf->A);
  *(intf->P) = P;      
    
  return err;
}

int solve_linear_system(void *m)
{
  int err = 0;
  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;
  
  intf->x = 0.0;
  itl::cyclic_iteration<double> iter(intf->b, 100, 1.e-11, 0.0, 5);
  cg(intf->AA, intf->xx, intf->bb, *(intf->P), iter);

  std :: cout << intf->xx <<"\n";
  return err;
}

void print_hello(void)
{
  printf("hello\n");
}


