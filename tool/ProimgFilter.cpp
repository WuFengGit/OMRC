// HogExtract.cpp : Defines the entry point for the console application.
//
#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"
#include <string.h>
using namespace std;
void find(int** bg, int bg_h, int bg_w, int** tar, int tar_h, int tar_w, int** curr_frame, int** samples, int len_sample, double*index)
{
    double** newpro = new double*[tar_h];
    for(int j=0;j<tar_h;j++)
        newpro[j] = new double[tar_w];
    for(int i=0;i<len_sample;i++)
    {
        for(int u=0;u<tar_w;u++)
            for(int v=0;v<tar_h;v++)
            {
                int x = samples[i][0]-1+u, y=samples[i][1]-1+v;
                if(bg[y][x]!=-1)
                {
                    int diff_curr_tar = 0, diff_curr_bg = 0;
                    for(int m=-2;m<=2;m++)
                        for(int n=-2;n<=2;n++)
                            if(v+m>=0&&v+m<tar_h&&u+n>=0&&u+n<tar_w)
                            {
                                if(bg[y+m][x+n]==-1)
                                    diff_curr_bg += abs(bg[y][x]-curr_frame[y+m][x+n]);
                                else
                                    diff_curr_bg += abs(bg[y+m][x+n]-curr_frame[y+m][x+n]);
                                
                                diff_curr_tar += abs(tar[v+m][u+n]-curr_frame[y+m][x+n]);
                            }
                    if(diff_curr_tar<=diff_curr_bg)
                        newpro[v][u] = 1;
                    else
                        newpro[v][u] = 0;
                    
                }else
                    newpro[v][u] = 1;
                
            }
        double sumpro = 0;
        for(int u=0;u<tar_w;u++)
            for(int v=0;v<tar_h;v++)
                sumpro += newpro[v][u];
        if(sumpro>tar_h*tar_w*0.5)
        {
            for(int j=0;j<tar_h;j++)
                delete newpro[j];
            delete newpro;
            index[0] = i+1;
            return;
        }  
    }
    for(int v=0;v<tar_h;v++)
        delete newpro[v];
    delete newpro;
    index[0] = 0;
    return;
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[])
        
{
    //update(int** bg, int bg_h, int bg_w, int** tar, int tar_h, int tar_w, int** curr_frame, int tar_y, int tar_x, int** pro)
    double* bg0, *tar0, *curr_frame0, *samples0;
    double* index;
    
    mwSize bg_h, bg_w, tar_h, tar_w, len_sample;
    
    bg0 = mxGetPr(prhs[0]);
    bg_h = mxGetM(prhs[0]);
    bg_w = mxGetN(prhs[0]);
    
    tar0 = mxGetPr(prhs[1]);
    tar_h = mxGetM(prhs[1]);
    tar_w = mxGetN(prhs[1]);
    
    curr_frame0 = mxGetPr(prhs[2]);
    samples0 = mxGetPr(prhs[3]);
    len_sample = mxGetM(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    index = mxGetPr(plhs[0]);
    
    int** bg = new int*[bg_h];
    int** curr_frame = new int*[bg_h];
    for(int i=0;i<bg_h;i++)
    {
        bg[i] = new int[bg_w];
        curr_frame[i] = new int[bg_w];
        for(int j=0;j<bg_w;j++)
        {
            bg[i][j] = bg0[i+bg_h*j];
            curr_frame[i][j] = curr_frame0[i+bg_h*j];
        }
    }
    
    int** tar = new int*[tar_h];
    for(int i=0;i<tar_h;i++)
    {
        tar[i] = new int[tar_w];
        for(int j=0;j<tar_w;j++)
        {
            tar[i][j] = tar0[i+tar_h*j];
        }
    }
    
    int** samples = new int*[len_sample];
    for(int i=0;i<len_sample;i++)
    {
        samples[i] = new int[4];
        for(int j=0;j<4;j++)
        {
            samples[i][j] = (int)samples0[i+len_sample*j];
        }
    }
 
    find(bg, bg_h, bg_w, tar, tar_h, tar_w, curr_frame, samples, len_sample, index);
    for(int i=0;i<bg_h;i++)
    {
        delete bg[i];
        delete curr_frame[i];
    }
    delete bg;
    bg = NULL;
    delete curr_frame;
    curr_frame = NULL;
    for(int i=0;i<tar_h;i++)
    {
        delete tar[i];
    }
    delete tar;
    tar = NULL;
    for(int i=0;i<len_sample;i++)
        delete samples[i];
    delete samples;
    samples = NULL;
    return;
    
}