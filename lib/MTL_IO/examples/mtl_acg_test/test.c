#include <stdio.h>
#include <stdlib.h>

#include "mtl_linear_solver_interface.h"

int main(int argc,char *argv[])
{
  int err = 0;
  void *intf = NULL;
  
  int elemno = 4;
  int nne = 4;
  int N = elemno*nne - 2*(elemno-1);
  
  printf("total number of degree freedom: %dx%d \n", N,N);
    
  double *b = (double *) malloc(sizeof(double)*N);  
  
  for(int I=0; I<N; I++)
    b[I] = 1.0;
  
  err += construct_solver(&intf);
  
  err += initialize_linear_system(intf, N);
  err += set_linear_system_b(intf, b, N); 
  
  
  for(int a=0; a<elemno; a++)
  {
    int id_0 = a*(nne - 2);
    int I[4] = {id_0, id_0+1, id_0+2, id_0+3};
    int J[4] = {id_0, id_0+1, id_0+2, id_0+3};

    int IN = 4;
    int JN = 4;
    double values[16] = {1.0,0.0,0.0,0.0,
                         0.0,1.0,0.0,0.0,
                         0.0,0.0,1.0,0.0,
                         0.0,0.0,0.0,1.0};
  
    err += update_linear_system_A_IJ(intf,I,IN,J,JN,values);
  }
                                          
  err += set_solver_pc(intf, 0);
  err += solve_linear_system(intf);
  
  err += destruct_solver(&intf);  
    
  printf("this is done\n");  
  free(b);  
  return 0;
}

