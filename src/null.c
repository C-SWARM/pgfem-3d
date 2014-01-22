#include "null.h"

void null_4d (double A[3][3][3][3])
{
  long i,j,k,l;
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      for (k=0;k<3;k++){
	for (l=0;l<3;l++){
	  A[i][j][k][l] = 0.0;
	}
      }
    }
  }  
}

void null_2d (double A[3][3])
{
  long i,j;
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      A[i][j] = 0.0;
    }
  }
}
