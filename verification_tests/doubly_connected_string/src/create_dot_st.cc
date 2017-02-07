#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;

int main(int argc,char *argv[])
{
  double tSteps, time1, time2, sum, dt;
  int t;
  tSteps = 3500;
	dt = 1e-8;

  time1 = 0;
  time2 = dt * tSteps;
  sum = time1;

  cout << "1.0e-6 5 3 1" << endl;
//how often to simulate
  cout << tSteps << endl;
  for (t = 0; t < tSteps + 1; t++) {
    cout << sum <<" ";
    sum = (t + 1)/tSteps * (time2 - time1);
  }
//How often to output
  int skip = 1;	
	cout << endl;
	cout << tSteps/skip << endl;
  for (t = 0; t < tSteps; t++) {
	  if(t%skip == 0) {
	    cout << t <<" ";
    } 	
  }

  cout << endl << "0" << endl; //no displacement steps

}
