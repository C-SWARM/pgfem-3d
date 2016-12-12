#include <stdio.h>
#include <stdlib.h> /* for rand() */
#include <unistd.h> /* for getpid() */
#include <time.h> /* for time() */
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>

using namespace std;

#define PI      3.1415926535897932384626433832795

int main(int argc,char *argv[])
{
        int i,j;
        int row = 0;
        int column;
        int nodes = 0;
        double d;
        double *positions;
        positions = new double[10000*3]; //assuming each chunk has less than 10k nodes
        string line;
        ifstream infile(argv[1]);
        while (getline(infile, line)) {
                istringstream iss(line);
                if(row == 0) {
                        if (iss >> d) {
                                nodes = d;
                        }
                }

                if((row > 7) && (row < (nodes + 8)) ) {
                        for (column = 0 ; column < 10; column++) {
                                if (iss >> d) {
                                        if ((column > 2) && (column < 6)) {
                                                positions[(row - 8) + (column - 3)*3] = d;
                                        }
                                }
                        }
                }
                row++;
        }

        double x;
        double L = 10; //length of bar
        cout<<"-1"<<endl;
        cout<<".5"<<endl;
        cout<<"1e-9"<<endl; //density = 1?
        for (i = 0; i < nodes; i++) {
                x = positions[i + 3*(0)];
                cout<<i<<" 0 "<<1e-1*cos(PI*x/L)<<" 0 0 0 0"<<endl;// half bar 
//		cout<<i<<" 0 "<<1e-1*sin(2*PI*x/L)<<" 0 0 0 0"<<endl;//full bar
//                cout<<i<<" 0 0 0 0 0 0"<<endl;// no initial displacement
        }


}


