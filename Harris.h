// Harris1.h: interface for the CHarris class.
//
//////////////////////////////////////////////////////////////////////

#ifndef HARRIS_H
#define HARRIS_H

#include "mex.h"
#include "utils.h"
#include <vector>
#include <math.h>

using namespace std;

typedef struct point{
  unsigned int x; 
  unsigned int y;
  double value; //cluster variance
}CPoint;

class CHarris  
{
private:
    int r0,c0;
    double **GrayValue;//存放灰度值，8位 GrayValue[r][c]表示倒数r行c列的灰度值
public:
    vector<CPoint> MaxpointArray;
public:
	void MakeGauss2(double sigma, double **pdKernel, int pnWindowSize);
	void Gaussian(double sigma, double **pUnchSmooth1,double **pUnchSmooth2,double **pUnchSmooth3, double **out1,double **out2, double **out3);
	void DirGrad2(double **pnGradX , double **pnGradY,double **pnGradXY);
	void Harris(double sigma,int thresh,int radius);
    void Convolve2D(double **image, double ** pnKernel, int pnWindowSize, double ** out);
public:
	CHarris(int row, int col, double **garyvalue);
	virtual ~CHarris();

};

#endif 
