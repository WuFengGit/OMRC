#include <cmath>
#include <ctime>
#include "mex.h"
using namespace std;

static void mask(double** proimg, int h, int w, int** rects, int N, double* bimask, double th)
{
    double** proimgIH = new double*[h];  
    for(int i=0;i<h;i++)
    {
        proimgIH[i] = new double[w];
        for(int j=0;j<w;j++)
            proimgIH[i][j] = 0;
    }
    proimgIH[0][0] = proimg[0][0];
    for(int i=1;i<h;i++)
        proimgIH[i][0] = proimgIH[i-1][0]+proimg[i][0];
    
    for(int j=1;j<w;j++)
    {
        proimgIH[0][j] = proimgIH[0][j-1]+proimg[0][j];
        for(int i=1;i<h;i++)
            proimgIH[i][j] = proimgIH[i-1][j]+proimgIH[i][j-1]-proimgIH[i-1][j-1]+proimg[i][j];
    }
    double* rects_avgpro = new double[N];
    for(int iter=0;iter<N;iter++)
    {
        int y2 = rects[iter][1]+rects[iter][3]-2;
        int x2 = rects[iter][0]+rects[iter][2]-2;
        int y1 = rects[iter][1]-2;
        int x1 = rects[iter][0]-2;
        if(x1>1&&y1>1)
            rects_avgpro[iter] = proimgIH[y2][x2] - proimgIH[y2][x1] -
                    proimgIH[y1][x2] + proimgIH[y1][x1];
        
        else if(x1>1)
            rects_avgpro[iter] = proimgIH[y2][x2] - proimgIH[y2][x1];
        
        else if(y1>1)
            rects_avgpro[iter] = proimgIH[y2][x2] - proimgIH[y1][x2];
        
        else
            rects_avgpro[iter] = proimgIH[y2][x2];
        
        rects_avgpro[iter] = rects_avgpro[iter]/rects[iter][2]/rects[iter][3];
    }
    
    int usec = 0;
    for(int i=0;i<N;i++)
        if(rects_avgpro[i]>=th)
        {
            bimask[i] = 1;
            usec ++;
        }
        else
            bimask[i] = 0;

    if(N<=25)
        for(int i=0;i<N;i++)
            bimask[i] = 1;
    else
        for(int i=0;i<25-usec;i++)
        {
            int maxidx = -1;
            for(int j=0;j<N;j++)
                if(bimask[j]==0&&(maxidx==-1||rects_avgpro[j]>rects_avgpro[maxidx]))
                    maxidx = j;
            bimask[maxidx] = 1;
        }
    
    delete rects_avgpro;
    rects_avgpro = NULL;
    for(int i=0;i<h;i++)
        delete proimgIH[i];
    delete proimgIH;
    proimgIH = NULL;
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    double*tproimg = mxGetPr(prhs[0]);
    double*trects = mxGetPr(prhs[1]);
    double*tth = mxGetPr(prhs[2]);
    
    size_t N, h, w;
    h = mxGetM(prhs[0]);
    w = mxGetN(prhs[0]);
    N = mxGetM(prhs[1]);
    double**proimg = new double*[h];
    for(int i=0;i<h;i++)
    {
        proimg[i] = new double[w];
        for(int j=0;j<w;j++)
            proimg[i][j] = tproimg[j*h+i];
    }
    int**rects = new int*[N];
    for(int i=0;i<N;i++)
    {
        rects[i] = new int[4];
        for(int j=0;j<4;j++)
            rects[i][j] = (int)trects[j*N+i];
    }
    double th = tth[0];
    plhs[0] = mxCreateDoubleMatrix(N,1, mxREAL);
    double* bimask = mxGetPr(plhs[0]);
    mask(proimg, h, w, rects, N, bimask, th);
    
    for(int i=0;i<h;i++)
        delete proimg[i];
    delete proimg;
    proimg = NULL;
    for(int i=0;i<N;i++)
        delete rects[i];
    delete rects;
    rects = NULL;
            
    return;
    
}