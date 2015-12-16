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
  itl::pc::identity<Distribute(CMatrix(double))> *P;

  char PC[1024];
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
//  MTL_SOLVER_INTF *temp = (MTL_SOLVER_INTF *) malloc(sizeof(MTL_SOLVER_INTF));
  MTL_SOLVER_INTF *temp = new MTL_SOLVER_INTF;
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

  delete temp;
  return err;
}

int initialize_linear_system(void *m, int N, int nne)
{ 
  int err = 0;  

  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;

  Vector(double) b(nne,0.0);
  CMatrix(double) A(nne,nne);
  A=0.0;
  DistributeV(Vector(double)) Pb(N,0.0);
  Distribute(CMatrix(double)) PA(N,N);
  PA=0.0;

  intf->A=A;
  intf->b=b;
  intf->x=b;

  intf->AA=PA;

  intf->bb=Pb;
  intf->xx=Pb;



  sprintf(intf->PC, "%s", "identity");
  return err;
}

int update_linear_system_A_IJ(void *m, int *I, int IN,
                                       int *J, int JN, int N, double *values, MPI_Comm comm_1)
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
    vI[a] = a;
    
  for(int a=0; a<JN; a++)
    vJ[a] = a;
  
  ins << element_matrix(temp, vI, vJ);
}
 // std :: cout << intf->A <<"\n";



  {
  mtl::matrix::inserter<Distribute(CMatrix(double)), mtl::operations::update_plus<double> > insA(intf->AA);
  mtl::vector::inserter< DistributeV(Vector(double)), mtl::operations::update_plus<double> > insV(intf->bb);

  for (int i=0;i<IN;i++){
	  int irow=I[i];
	  insV[irow] << intf->b[i];
     for (int j=0;j<JN;j++) {
      int icol=J[i];
       insA[irow][icol] << intf->A[i][j];
     }
    }
  }


    mtl::par::sout << "The Matrix AA is\n" << intf->AA << "\n";
    mtl::par::sout << "The Vector bb is\n" << intf->bb << "\n";

  return err;
}



int set_linear_system_b(void *m, double *values, int N)
{   
  int err = 0;  


  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;
  boost::mpi::communicator comm(communicator(intf->AA));

  double *b = (intf->b).address_data();
        
  for(int a=0; a<N; a++)
   b[a] = values[a];
  
  return err;
}

int set_solver_pc(void *m, int type)
{   
  int err = 0;

  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;  
  itl::pc::identity<Distribute(CMatrix(double))> P(intf->AA);
 // itl::pc::ilu_0<CMatrix(double)>    P(intf->A);

  *(intf->P) = P;
    
  return err;
}

int solve_linear_system(void *m,int N)
{
  int err = 0;
  MTL_SOLVER_INTF *intf = (MTL_SOLVER_INTF *) m;
  
   intf->xx = 0.0;

    itl::cyclic_iteration<double> iter(intf->bb, 500, 1.e-6, 0.0, 5);

    itl::pc::identity<Distribute(CMatrix(double))> PC(intf->AA);

   if (intf->PC=="ilu_0") {
    itl::pc::ilu_0<Distribute(CMatrix(double))> PC(intf->AA);
    }
    else if (intf->PC=="diagonal") {
    itl::pc::diagonal<Distribute(CMatrix(double))> PC(intf->AA);
    }

      bicgstab_2(intf->AA, intf->xx, intf->bb, PC, iter);

  return err;
}

void print_hello(void)
{
  printf("hello\n");
}


