// HogExtract.cpp : Defines the entry point for the console application.
//
#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"
#include <string.h>
#include <time.h>
const int MAX_MOV = 25;
const int R = 2;
const int HIGH = 30;
const int LOW = 10;
const int DILATE_ERODE = 3;

void gaussianFilter(double** mat, int mat_h, int mat_w, double sigma, double**& output);
void imdilate_erode(int** pro, int h, int w);

double findTH(int** curr_frame_crop, int** tar, int** th_pro,double th0, int tar_h, int tar_w)
{
    int bincount_w = 2, bincount_h = 2;
    int binsize_w = tar_w/bincount_w, binsize_h = tar_h/bincount_h;
    double min_sigma = 255;
    double min_avg = 255;
    
    double** diff = new double*[tar_h];
    for(int v=0;v<tar_h;v++)
    {
        diff[v] = new double[tar_w];
        for(int u=0;u<tar_w;u++)
            diff[v][u] = fabs(curr_frame_crop[v][u]-tar[v][u]);
    }
    
    for(int v = 0;v<bincount_h;v++)
        for(int u=0;u<bincount_w;u++)
        {
            double avg = 0, sigma = 0;
            int count = 0;
            for(int m=0;m<binsize_h;m++)
                for(int n=0;n<binsize_w;n++)
                {
                    avg += diff[v*binsize_h+m][u*binsize_w+n];
                    count ++;
                }                    
            avg /= count;
            for(int m=0;m<binsize_h;m++)
                for(int n=0;n<binsize_w;n++)
                {
                    sigma += pow(diff[v*binsize_h+m][u*binsize_w+n]-avg,2);
                }
            sigma = sqrt(sigma/count);
            if(min_sigma+min_avg>sigma+avg)
            {
                min_sigma = sigma;
                min_avg = avg;
            }
        }
    
    //min_sigma *= 2; 
    if(th0>0)
        th0 = (min_avg+min_sigma)*0.2+th0*0.8;
    else
        th0 = min_avg+min_sigma;
        

    // printf("min_sigma %f\n", min_sigma);
    for(int v=0;v<tar_h;v++)
        for(int u=0;u<tar_w;u++)
            if(diff[v][u] <= th0)
                th_pro[v][u] = 1;
            else
                th_pro[v][u] = 0;
    
    for(int v=0;v<tar_h;v++)
        delete diff[v];
    delete diff;
    diff = NULL;
    return th0;
}

void gaussianFilter(double** mat, int mat_h, int mat_w, double sigma, double**& output)
{
    double** filter = new double*[2*R+1];
    double filtersum = 0;
    for(int i=0;i<2*R+1;i++)
    {
        filter[i] = new double[2*R+1];
        for(int j=0;j<2*R+1;j++)
        {
            filter[i][j] = exp(-(pow(i-R,2)+pow(j-R,2))/2/pow(sigma,2));
            filtersum += filter[i][j];
        }
    }
    for(int i=0;i<2*R+1;i++)
        for(int j=0;j<2*R+1;j++)
            filter[i][j] /= filtersum;
    
    double** tmpmat = new double*[mat_h+2*R];
    for(int i=0;i<mat_h+2*R;i++)
    {
        tmpmat[i] = new double[mat_w+2*R];
        for(int j=0;j<mat_w+2*R;j++)
        {
            int m = i<=R?0:i-R-1;
            m = m>=mat_h?mat_h-1:m;
            int n = j<=R?0:i-R-1;
            n = n>=mat_w?mat_w-1:n;
            tmpmat[i][j] = mat[m][n];            
        }
    }
    
    for(int i=0;i<mat_h;i++)
        for(int j=0;j<mat_w;j++)
        {
            for(int m=0;m<2*R+1;m++)
                for(int n=0;n<2*R+1;n++)
                    output[i][j] += tmpmat[i+m][j+n]*filter[m][n];
        }
    
    // delete
    for(int i=0;i<2*R+1;i++)
        delete filter[i];
    delete filter;
    filter = NULL;
    for(int i=0;i<mat_h+2*R;i++)
        delete tmpmat[i];
    delete tmpmat;
    tmpmat = NULL;            
}

double update(int** bg, int bg_h, int bg_w, int** tar, int tar_h, int tar_w, int** curr_frame, int tar_y, int tar_x, int** pro, double th0, double usePro)
{
    //clock_t start;
    //start = clock();
    //printf("%d\t", clock()-start);
    float**dist = new float*[MAX_MOV*2+1];
    for(int i=0;i<MAX_MOV*2+1;i++)
        dist[i] = new float[MAX_MOV*2+1];
    
   
    
    // full image
    /*for(int v=-MAX_MOV;v<=MAX_MOV;v++)
        for(int u=-MAX_MOV;u<=MAX_MOV;u++)
        {
            int t_count = 0;
            float t_dist = 0;
            for(int i=0;i<bg_h;i++)
                for(int j=0;j<bg_w;j++)
                    if(bg[i][j]!=-1&& i-v>=0 && i-v<bg_h&& j-u>=0&&j-u<bg_w&&
                     !(i-v>=tar_y&&i-v<tar_y+tar_h&&j-u>=tar_x&&j-u<tar_x+tar_w))
                    {
                        t_count ++;
                        t_dist += abs(bg[i][j]-curr_frame[i-v][j-u]);
                    }
            dist[v+MAX_MOV][u+MAX_MOV] = t_dist/t_count;
        }*/
    
    
    for(int v=-MAX_MOV;v<=MAX_MOV;v++)
        for(int u=-MAX_MOV;u<=MAX_MOV;u++)
        {
            int t_count = 0;
            float t_dist = 0;
            int istart = tar_y-2.5*tar_h+v;
            istart=istart>0?istart:0;
            int iend = tar_y+3.5*tar_h+v;      
            iend=iend<tar_h?iend:tar_h;    
            int jstart = tar_x-2.5*tar_w+u;
            jstart=jstart>0?jstart:0;
            int jend = tar_x+3.5*tar_w+u; 
            jend=jend<tar_w?jend:tar_w;
            
            for(int i=istart;i<iend;i++)
                for(int j=jstart;j<jend;j++)                
                    if(bg[i][j]!=-1&& i-v>=0 && i-v<bg_h&& j-u>=0&&j-u<bg_w&&
                     !(i-v>=tar_y&&i-v<tar_y+tar_h&&j-u>=tar_x&&j-u<tar_x+tar_w))
                    {
                        t_count ++;
                        t_dist += abs(bg[i][j]-curr_frame[i-v][j-u]);
                    }
            dist[v+MAX_MOV][u+MAX_MOV] = t_dist/t_count;
        }
    
    // img1 & img2
    /*
    int **c_img1 = new int*[bg_h+MAX_MOV*2+1];
    int **c_img2 = new int*[bg_h+MAX_MOV*2+1];
    for(int i=0;i<bg_h+MAX_MOV*2+1;i++)
    {
        c_img1[i] = new int[bg_w+MAX_MOV*2+1];
        c_img2[i] = new int[bg_w+MAX_MOV*2+1];
        for(int j=0;j<bg_w+MAX_MOV*2+1;j++)
            c_img1[i][j] = -1;
    }

    for(int i=0;i<bg_h;i++)
        for(int j=0;j<bg_w;j++)
            c_img1[MAX_MOV+i][MAX_MOV+j] = bg[i][j];
    for(int v=-MAX_MOV;v<=MAX_MOV;v++)
        for(int u=-MAX_MOV;u<=MAX_MOV;u++)
        {
            for(int i=0;i<bg_h+MAX_MOV*2+1;i++)
                for(int j=0;j<bg_w+MAX_MOV*2+1;j++)
                    c_img2[i][j] = -1;
            for(int i=0;i<bg_h;i++)
                for(int j=0;j<bg_w;j++)
                    c_img2[i+v+MAX_MOV][j+u+MAX_MOV] = curr_frame[i][j];
            for(int i=0;i<tar_h;i++)
                for(int j=0;j<tar_w;j++)
                    c_img2[i+v+MAX_MOV+tar_y-1][j+u+MAX_MOV+tar_x-1] = -1;
            int t_count = 0;
            float t_dist = 0;
            for(int i=0;i<bg_h+2*MAX_MOV+1;i++)
                for(int j=0;j<bg_w+2*MAX_MOV+1;j++)
                    if(c_img1[i][j]!=-1&&c_img2[i][j]!=-1)
                    {
                        t_count ++;
                        t_dist += abs(c_img1[i][j]-c_img2[i][j]);
                    }
            dist[v+MAX_MOV][u+MAX_MOV] = t_dist/t_count;
        }*/
    int min_dist = 255, min_mov_y=0, min_mov_x = 0;
    for(int i=0;i<2*MAX_MOV+1;i++)
        for(int j=0;j<2*MAX_MOV+1;j++)
            if(min_dist>dist[i][j])
            {
                min_dist = dist[i][j];
                min_mov_y = -(i-MAX_MOV);
                min_mov_x = -(j-MAX_MOV);
            }

    for(int i=0;i<MAX_MOV*2+1;i++)
        delete dist[i];
    delete dist;
    dist = NULL;
    /*for(int i=0;i<bg_h+MAX_MOV*2+1;i++)
    {
        delete c_img1[i];
        delete c_img2[i];
    }
    delete c_img1;
    delete c_img2;
    c_img1 = NULL;
    c_img2 = NULL;*/

    // update bgmodel
    int ** newbg = new int*[bg_h];
    for(int i=0;i<bg_h;i++)
        newbg[i] = new int[bg_w];
    for(int v=0;v<bg_h;v++)
        for(int u=0;u<bg_w;u++)
            if(v-min_mov_y>=0&&v-min_mov_y<bg_h&&u-min_mov_x>=0&&u-min_mov_x<bg_w && bg[v-min_mov_y][u-min_mov_x]!=-1)
            {
                if(v>=tar_y-1&&v<tar_y-1+tar_h&&u>=tar_x-1&&u<tar_x-1+tar_w)
                    newbg[v][u] = bg[v-min_mov_y][u-min_mov_x];
                else
                    newbg[v][u] = bg[v-min_mov_y][u-min_mov_x]*0.85+curr_frame[v][u]*0.15;
            }else{
                if(v>=tar_y-1&&v<tar_y-1+tar_h&&u>=tar_x-1&&u<tar_x-1+tar_w)
                    newbg[v][u] = -1;
                else
                    newbg[v][u] = curr_frame[v][u];
            }
    for(int v=0;v<bg_h;v++)
        for(int u=0;u<bg_w;u++)
            bg[v][u] = newbg[v][u];

    for(int i=0;i<bg_h;i++)
        delete newbg[i];
    delete newbg;
    newbg = NULL;
    
    // update pro
    int** newpro = new int*[tar_h];
    for(int i=0;i<tar_h;i++)
        newpro[i] = new int[tar_w];
    
    for(int u=0;u<tar_w;u++)
        for(int v=0;v<tar_h;v++)
        {
            int x = tar_x-1+u, y=tar_y-1+v;
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
            }
            else
               newpro[v][u] = 1; 
        }
    int ** curr_frame_crop = new int*[tar_h];
    int ** th_pro = new int*[tar_h];
    for(int v=0;v<tar_h;v++)
    {
        th_pro[v] = new int[tar_w];
        curr_frame_crop[v] = new int[tar_w];
        for(int u=0;u<tar_w;u++)
            curr_frame_crop[v][u] = curr_frame[tar_y-1+v][tar_x-1+u];
    }
    
    double th1 = findTH(curr_frame_crop, tar, th_pro, th0, tar_h, tar_w);
    for(int u=0;u<tar_w;u++)
        for(int v=0;v<tar_h;v++)
            if(bg[tar_y-1+v][tar_x-1+u]==-1)
                newpro[v][u] = th_pro[v][u];
    for(int u=0;u<tar_w;u++)
        for(int v=0;v<tar_h;v++)
            pro[v][u] = newpro[v][u];
    
    for(int v=0;v<tar_h;v++)
        delete newpro[v];
    delete newpro;
    //imdilate_erode(pro, tar_h, tar_w);
    if(usePro<1)
        for(int u=0;u<tar_w;u++)
        for(int v=0;v<tar_h;v++)
            pro[v][u] = 1;
    
    // update tar
    for(int u=0;u<tar_w;u++)
        for(int v=0;v<tar_h;v++)
            if(pro[v][u] == 1)
                tar[v][u] = tar[v][u]*0.85+curr_frame[v+tar_y-1][u+tar_x-1]*0.15;
    return th1;
}
void imdilate_erode(int** pro, int h, int w)
{
    int** tmp = new int*[h];
    for(int i=0;i<h;i++)
    {
        tmp[i] = new int[w];
        for(int j=0;j<w;j++)
            tmp[i][j] = pro[i][j];
    }
    for(int i=0;i<h;i++)
        for(int j=0;j<w;j++)
            if(pro[i][j]==1)
            {
                for(int m=-DILATE_ERODE/2;m<=DILATE_ERODE/2;m++)
                    for(int n=-DILATE_ERODE/2;n<=DILATE_ERODE/2;n++)
                        if(m+i>=0&&n+j>=0&&m+i<h&&n+j<w)
                            tmp[i+m][j+n] = 1;
            }
    for(int i=0;i<h;i++)
        for(int j=0;j<w;j++)
            pro[i][j] = tmp[i][j];
    for(int i=0;i<h;i++)
        for(int j=0;j<w;j++)
        {
            if(tmp[i][j]==0)
            {
                for(int m=-DILATE_ERODE/2;m<=DILATE_ERODE/2;m++)
                    for(int n=-DILATE_ERODE/2;n<=DILATE_ERODE/2;n++)
                        if(m+i>=0&&n+j>=0&&m+i<h&&n+j<w)
                            pro[i+m][j+n] = 0 ;
            }
        }   
    
    for(int i=0;i<h;i++)
        delete tmp[i];
    delete tmp;
    tmp = NULL;
}



void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[])
        
{
    double* bg0, *tar0, *curr_frame0, *pro0, *tar_y0, *tar_x0, *th0;
    double* bg1, *tar1, *pro1, *th1;
    
    mwSize bg_h, bg_w, tar_h, tar_w;
    
    bg0 = mxGetPr(prhs[0]);
    bg_h = mxGetM(prhs[0]);
    bg_w = mxGetN(prhs[0]);
    
    tar0 = mxGetPr(prhs[1]);
    tar_h = mxGetM(prhs[1]);
    tar_w = mxGetN(prhs[1]);
    
    curr_frame0 = mxGetPr(prhs[2]);
    pro0 = mxGetPr(prhs[3]);
    tar_y0 = mxGetPr(prhs[4]);
    tar_x0 = mxGetPr(prhs[5]);
    th0 = mxGetPr(prhs[6]);
    
    plhs[0] = mxCreateDoubleMatrix(bg_h, bg_w, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(tar_h, tar_w, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(tar_h, tar_w, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    bg1 = mxGetPr(plhs[0]);
    tar1 = mxGetPr(plhs[1]);
    pro1 = mxGetPr(plhs[2]);
    th1 = mxGetPr(plhs[3]);
    
    //update(int** bg, int bg_h, int bg_w, int** tar, int tar_h, int tar_h, int** curr_frame, int tar_y, int tar_x, int** pro)

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
    int** pro = new int*[tar_h];
    for(int i=0;i<tar_h;i++)
    {
        tar[i] = new int[tar_w];
        pro[i] = new int[tar_w];
        for(int j=0;j<tar_w;j++)
        {
            tar[i][j] = tar0[i+tar_h*j];
            pro[i][j] = pro0[i+tar_h*j];
        }
    }
    
    int tar_y = (int)(tar_y0[0]);
    int tar_x = (int)tar_x0[0];
    th1[0] = update(bg, bg_h, bg_w, tar, tar_h, tar_w, curr_frame, tar_y, tar_x, pro, th0[0], th0[1]);
    
    for(int i=0;i<bg_h;i++)
        for(int j=0;j<bg_w;j++)
            bg1[i+bg_h*j] = bg[i][j];
    
    for(int i=0;i<tar_h;i++)
        for(int j=0;j<tar_w;j++)
        {
            tar1[i+tar_h*j] = tar[i][j];
            pro1[i+tar_h*j] = pro[i][j];
        }
    
    
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
        delete pro[i];
    }
    delete tar;
    tar = NULL;
    delete pro;
    pro = NULL;
    return;
    
}