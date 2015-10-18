#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;

double **calloc_mat(int dimX, int dimY){
    double **m = (double **)calloc(dimX, sizeof(double*));
    double *p = (double *)calloc(dimX*dimY, sizeof(double));
    
    for(int i=0; i <dimX;i++){
    m[i] = &p[i*dimY];

    }
   return m;
}

void free_mat(double **m){
  free(m[0]);
  free(m);
}

void write_mat_file(int dimX, int dimY, double **mat, const char *name)
{
    ofstream outdata; 
    outdata.open(name); 
    if( !outdata ) 
    { // file couldn't be opened
      exit(1);
    }
    for(int i=0; i<dimX; i++)
    {
        for(int j=0; j<dimY; j++)
        {
            outdata<<mat[i][j];
            if(j!=dimY-1)
                outdata<<" ";
        }

        outdata<<endl;
    }
    outdata.close();
}

#endif 