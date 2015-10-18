#include "Harris.h"
#include <iostream>

bool debug=true;
using namespace std;

CHarris::CHarris(int row, int col, double **garyvalue)
{
    r0=row;
    c0=col;
    GrayValue=garyvalue;
}
CHarris::~CHarris()
{
    MaxpointArray.clear();
}


void CHarris::Convolve2D(double **image, double ** pnKernel, int pnWindowSize, double **out)
{
	int kCenterX=pnWindowSize/2;
    int kCenterY=pnWindowSize/2;

    int rowIndex,colIndex;
	int mm=0,nn=0;

    for (int i=0;i<r0;i++)
        for (int j=0;j<c0;j++)
        {
            for (int m=0;m<pnWindowSize;m++)
            {
				mm=pnWindowSize-1-m;
                for (int n=0;n<pnWindowSize;n++)
                {
                    // index of input signal,used for cheking bundary
					nn=pnWindowSize-1-n;
                    rowIndex=i+(m-kCenterY);
                    colIndex=j+(n-kCenterX);
                    //ignore input samples which are out of bound
                    if (rowIndex>=0&&rowIndex<r0&&colIndex>=0&&colIndex<c0)
                    {
                       out[i][j]+=image[rowIndex][colIndex]* pnKernel[mm][nn];
                    }
                }

            }
        }
}
void CHarris::Harris(double sigma,int thresh,int radius)
{

	double ** pnGradX ; 
	double ** pnGradY ;
	double ** pnGradXY ;
	double ** pMaxValue;
	int y ,i,j;
	int x ; int m_maxvalue;
	CPoint maxpoint;//角点

	MaxpointArray.clear();
	pnGradX      = calloc_mat(r0,c0);
	pnGradY      = calloc_mat(r0,c0);
	pnGradXY     = calloc_mat(r0,c0);
	pMaxValue    = calloc_mat(r0,c0);
	// 计算方向导数
	DirGrad2(pnGradX, pnGradY,pnGradXY) ;
	
	// 对gx gy gxgy 进行高斯滤波
	double **out1,**out2,**out3;

    out1 = calloc_mat(r0,c0);
    out2 = calloc_mat(r0,c0);
    out3 = calloc_mat(r0,c0);
	Gaussian(sigma, pnGradX,pnGradY,pnGradXY,out1,out2,out3) ;
	//求R=(AB-C^2)-k(A+B)^2 
    double maxValue=0.0;
    
	for(y=0; y<(int)r0; y++)
		for(x=0; x<(int)c0; x++)
		{
            pMaxValue[y][x]=(out1[y][x]*out2[y][x]-out3[y][x]*out3[y][x])\
                            -0.04*(out1[y][x]+out2[y][x])*(out1[y][x]+out2[y][x]);
            if(pMaxValue[y][x]>maxValue)
                maxValue=pMaxValue[y][x];   // Find the maxmium of pMaxValue
		}

	// Make every value between 0-255
    maxValue=1000/maxValue;
    for(y=0; y<(int)r0; y++)
		for(x=0; x<(int)c0; x++)
		{
            pMaxValue[y][x]=maxValue*pMaxValue[y][x];
		}
	
    //write_mat_file(r0,c0,pMaxValue,"pmaxvalue.txt");

    
	for(int y=0; y<(int)r0; y+=radius)
	{
		for(int x=0; x<(int)c0; x+=radius)
		{
			m_maxvalue=-1000000;
            if(y>=int(radius) && y<r0-int(radius) && x>=int(radius) && x<c0-int(radius))
            {
                for(int m=-int(radius); m<=int(radius); m++)
                    for(int n=-int(radius); n<=int(radius); n++)
                {
                    if(pMaxValue[y+m][x+n]>m_maxvalue)
                        m_maxvalue=(int)pMaxValue[y+m][x+n];
                        maxpoint.x=y+m;
                        maxpoint.y=x+n;
                }
                for(int m=-int(radius); m<=int(radius); m++)
                    for(int n=-int(radius); n<=int(radius); n++)
                {
                    if(pMaxValue[y+m][x+n]!=m_maxvalue)
                        pMaxValue[y+m][x+n]=0;
                }
                if(m_maxvalue>thresh) 
                {
                    MaxpointArray.push_back(maxpoint);//存入角点坐标
                }

            }
		}
	}
    write_mat_file(r0,c0,pMaxValue,"pmaxvalue.txt");
	// 释放内存
	free_mat(pnGradX);
	free_mat(pnGradY);
	free_mat(pnGradXY);

    free_mat(out1);
    free_mat(out2);
    free_mat(out3);

}

void CHarris::DirGrad2(double **pnGradX , double **pnGradY, double **pnGradXY)
{
    double **Ix=calloc_mat(3,3);
    double **Iy=calloc_mat(3,3);
	
	//Ix=[-1 0 1;-1 0 1;-1 0 1] Iy=Ix'
    for(int i=0; i<3; i++)
	{
        for(int j=0; j<3; j++)
		{
			if(j==0)
				Ix[i][j]=-1;
			if(i==0)
				Iy[i][j]=-1;
			if(j==1)
				Ix[i][j]=0;
			if(i==1)
				Iy[i][j]=0;
			if(j==2)
				Ix[i][j]=1;
			if(i==2)
				Iy[i][j]=1;
		}
	}

	Convolve2D(GrayValue,Ix,3,pnGradX);
	Convolve2D(GrayValue,Iy,3,pnGradY);
    
	for(int y=0; y<(int)r0; y++)
		for(int x=0; x<(int)c0; x++)
			pnGradXY[y][x]=pnGradX[y][x]*pnGradY[y][x];
	for(int y=0; y<(int)r0; y++)
		for(int x=0; x<(int)c0; x++)
		{
			pnGradX[y][x]=pnGradX[y][x]*pnGradX[y][x];
			pnGradY[y][x]=pnGradY[y][x]*pnGradY[y][x];	
		}
        
}

void CHarris::Gaussian(double sigma, double **pUnchSmooth1,double **pUnchSmooth2,double **pUnchSmooth3, double **out1, double **out2, double **out3)
{
    int nWinSize=(int)floor(6*sigma);
	if(nWinSize%2==0)
		nWinSize++;		//Make nWindSize Odd!
    
    double ** Gtemplate=calloc_mat(nWinSize,nWinSize);
    
    MakeGauss2(sigma,Gtemplate,nWinSize);
    Convolve2D(pUnchSmooth1, Gtemplate, nWinSize,out1);
    Convolve2D(pUnchSmooth2, Gtemplate, nWinSize,out2);
    Convolve2D(pUnchSmooth3, Gtemplate, nWinSize,out3);
}

void CHarris::MakeGauss2(double sigma, double **pdKernel, int pnWindowSize)
{   //产生一个二维高斯模板，pnWindowsSize为模板的长和宽
	int i,j ;
	
	// 数组的中心点
	double nCenter;
	// 数组的某一点到中心点的距离
	double  dDis,dDiy  ; 

	//double PI = 3.14159;
	// 中间变量
	double  dValue; 
	double  dSum  ;
	dSum = 0 ; 
	// 中心
	nCenter = ((double)pnWindowSize-1) / 2;
	for(i=0; i<(pnWindowSize); i++)
	{
		for(j=0; j<(pnWindowSize); j++)
		{
			dDis = (double)(i - nCenter);
			dDiy = (double)(j - nCenter);
			dValue=exp(-(dDis*dDis+dDiy*dDiy)/(2*sigma*sigma));
			pdKernel[i][j]=dValue;
			dSum += dValue;
		}
	}
	for(i=0; i<(pnWindowSize); i++)//归一化
	{
		for(j=0; j<(pnWindowSize); j++)
		{
			pdKernel[i][j]/= dSum;
            pdKernel[i][j]=floor(pdKernel[i][j] * 10000.000f + 0.5) / 10000.000f; //保留小数点后面4位
		}
	}
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	if(nrhs!=4)
		mexErrMsgTxt("Wrong number of inputs");
	if(nlhs!=2)
		mexErrMsgTxt("Wrong number of outputs");
	double *img=(double *)mxGetPr(prhs[0]);
	double sigma=mxGetScalar(prhs[1]);
	int thresh=(int)mxGetScalar(prhs[2]);
	int radius=(int)mxGetScalar(prhs[3]);

	int row=mxGetM(prhs[0]);	//Row of image
	int col=mxGetN(prhs[0]);	//colum of image

	double **image=calloc_mat(row,col);

	for(int i=0; i<row; i++)
		for(int j=0; j<col; j++)
			image[i][j]=*(img+j*row+i);

	CHarris harris(row,col,image);
	harris.Harris(sigma,thresh,radius);

	plhs[0]=mxCreateDoubleMatrix(1,harris.MaxpointArray.size(),mxREAL);
	plhs[1]=mxCreateDoubleMatrix(1,harris.MaxpointArray.size(),mxREAL);

	double *x=(double *)mxGetData(plhs[0]);
	double *y=(double *)mxGetData(plhs[1]);
	int i=0;
	for( vector<CPoint>::iterator it=harris.MaxpointArray.begin(); it!=harris.MaxpointArray.end(); it++)
	{
		*(x+i)=(double)it->x;
		*(y+i)=(double)it->y;
		i++;
	}
	free_mat(image);	
}

