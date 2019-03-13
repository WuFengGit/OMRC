// HogExtract.cpp : Defines the entry point for the console application.
//
#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"
#include <string.h>

void maxminArray(int * arr, int len, int& maxV, int& minV)
{
	maxV = minV = arr[0];
	for (int i = 1; i < len; i++)
	{
		if (maxV < arr[i])
			maxV = arr[i];
		if (minV > arr[i])
			minV = arr[i];
	}
}

double*** IntegralImg(int ** im, int rows, int cols)
{
	double** Ix = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		Ix[i] = new double[cols];
		Ix[i][0] = Ix[i][cols-1] = 0;
		for (int j = 0; j < cols-2; j++)
			Ix[i][j+1] = im[i][j] - im[i][j +2];
	}
	double** Iy = new double*[rows];
	for (int i = 0; i < rows; i++)
        Iy[i] = new double[cols];
    for (int j = 0; j < cols; j++)
		Iy[0][j] = Iy[rows-1][j] = 0;
	for (int i = 0; i < rows-2; i++)
		for (int j = 0; j < cols; j++)
			Iy[i+1][j] = im[i][j] - im[i+2][j];
	
    double** angle = new double*[rows];
	double** magnitude = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		angle[i] = new double[cols];
		magnitude[i] = new double[cols];
		for (int j = 0; j < cols; j++)
		{
            magnitude[i][j] = sqrt(Ix[i][j]*Ix[i][j] + Iy[i][j]*Iy[i][j]);
            if (Iy[i][j]==0)
                angle[i][j] = 0;
            else
                angle[i][j] = atan(Ix[i][j] / Iy[i][j])/(M_PI/2)*90 + 90; 			
		}
	}
    
	for (int i = 0; i < rows; i++)
        delete Ix[i];
    delete Ix;
    Ix = NULL;
    for (int i = 0; i < rows; i++)
        delete Iy[i];
    delete Iy;
    Iy = NULL;

	double***iH = new double**[rows];
	for (int i = 0; i < rows; i++)
	{
		iH[i] = new double*[cols];
		for (int j = 0; j < cols; j++)
		{
			iH[i][j] = new double[9];
			for (int k = 0; k < 9; k++)
				iH[i][j][k] = 0;

			double alpha = angle[i][j];

			if (alpha > 10 && alpha <= 30)
			{
				iH[i][j][0] += magnitude[i][j] * (30 - alpha) / 20;
				iH[i][j][1] += magnitude[i][j]*(alpha - 10) / 20;
			}
			else if(alpha>30 && alpha <= 50)
			{
				iH[i][j][1] += magnitude[i][j] *(50 - alpha) / 20;
				iH[i][j][2] += magnitude[i][j] *(alpha - 30) / 20;
			}
			else if(alpha > 50 && alpha <= 70)
			{
				iH[i][j][2] += magnitude[i][j] *(70 - alpha) / 20;
				iH[i][j][3] += magnitude[i][j] *(alpha - 50) / 20;
			}
			else if(alpha > 70 && alpha <= 90)
			{
				iH[i][j][3] += magnitude[i][j] *(90 - alpha) / 20;
				iH[i][j][4] += magnitude[i][j] *(alpha - 70) / 20;
			}
			else if(alpha > 90 && alpha <= 110)
			{
				iH[i][j][4] += magnitude[i][j] *(110 - alpha) / 20;
				iH[i][j][5] += magnitude[i][j] *(alpha - 90) / 20;
			}
			else if (alpha>110 && alpha <= 130)
			{
				iH[i][j][5] += magnitude[i][j] *(130 - alpha) / 20;
				iH[i][j][6] += magnitude[i][j] *(alpha - 110) / 20;
			}
			else if (alpha>130 && alpha <= 150)
			{
				iH[i][j][6] += magnitude[i][j] *(150 - alpha) / 20;
				iH[i][j][7] += magnitude[i][j] *(alpha - 130) / 20;
			}
			else if (alpha>150 && alpha <= 170)
			{
				iH[i][j][7] += magnitude[i][j] *(170 - alpha) / 20;
				iH[i][j][8] += magnitude[i][j] *(alpha - 150) / 20;
			}
			else if (alpha >= 0 && alpha <= 10)
			{
				iH[i][j][0] += magnitude[i][j] *(alpha + 10) / 20;
				iH[i][j][8] += magnitude[i][j] *(10 - alpha) / 20;
			}
			else if (alpha>170 && alpha <= 180)
			{
				iH[i][j][8] += magnitude[i][j] *(190 - alpha) / 20;
				iH[i][j][0] += magnitude[i][j] *(alpha - 170) / 20;
			}
		}
	}
    
    for (int i = 0; i < rows; i++)
        delete angle[i];
    delete angle;
    angle = NULL;
    for (int i = 0; i < rows; i++)
        delete magnitude[i];
    delete magnitude;
    magnitude = NULL;

	for (int j = 1; j < cols; j++)
		for (int k = 0; k < 9;k++)
			iH[0][j][k] = iH[0][j][k] + iH[0][j-1][k];
	for (int i = 1; i < rows; i++)
		for (int k = 0; k < 9; k++)
			iH[i][0][k] = iH[i][0][k] + iH[i - 1][0][k];
    for (int i = 1; i < rows; i++)
		for (int j = 1; j < cols; j++)
			for (int k = 0; k < 9; k++)
				iH[i][j][k] = iH[i][j][k] + iH[i][j - 1][k] + iH[i - 1][j][k] - iH[i - 1][j - 1][k];    
    return iH;      
}

static void getFtrValHOG(int img0[], int w, int h, int samples0[], int sampleLen, int ftrpos0[], int ftrLen, double* hog0)
{
    int** img = new int*[h];
	for (int i = 0; i < h; i++)
	{
		img[i] = new int[w];
		memcpy(img[i], img0+i*w, w * sizeof(int));
	}

	int* samplex = samples0;
	int* sampley = samples0 + sampleLen;

	int* ftrposx = ftrpos0;
	int* ftrposy = ftrpos0 + ftrLen;
	int* ftrposw = ftrpos0 + ftrLen*2;
	int* ftrposh = ftrpos0 + ftrLen*3;

	// crop img
	int maxSx, minSx, maxSy, minSy;
	maxminArray(samplex, sampleLen, maxSx, minSx);
	maxminArray(sampley, sampleLen, maxSy, minSy);
	int maxFrtposx, maxFrtposy, maxFrtposxw, maxFrtposyh, tmpV, minFrtposx, minFrtposy;
	maxminArray(ftrposx, ftrLen, maxFrtposx, minFrtposx);
	maxminArray(ftrposy, ftrLen, maxFrtposy, minFrtposy);
	int* tmpArr = new int[ftrLen];
	for (int i = 0; i < ftrLen; i++)
		tmpArr[i] = ftrposx[i] + ftrposw[i]-1;
	maxminArray(tmpArr, ftrLen, maxFrtposxw, tmpV);
	for (int i = 0; i < ftrLen; i++)
		tmpArr[i] = ftrposy[i] + ftrposh[i]-1;
	maxminArray(tmpArr, ftrLen, maxFrtposyh, tmpV);
	int imgminx = minSx + minFrtposx - 2;
	int imgminy = minSy + minFrtposy - 2;
	int imgmaxx = maxSx + maxFrtposxw - 2;
	int imgmaxy = maxSy + maxFrtposyh - 2;
    delete tmpArr;
    tmpArr = NULL;
    
	int **newimg = new int*[imgmaxy - imgminy + 1];
	for (int i = 0; i < imgmaxy - imgminy + 1; i++)
	{
		newimg[i] = new int[imgmaxx - imgminx + 1];
		memcpy(newimg[i], img[i+imgminy] + imgminx, (imgmaxx - imgminx + 1) * sizeof(int));
	}
    for (int i = 0; i < h; i++)
        delete img[i];
    delete img;
    img = NULL;
	// intergeral image
	double*** iH = IntegralImg(newimg, imgmaxy - imgminy + 1, imgmaxx - imgminx + 1); 
    for(int i=0;i<imgmaxy - imgminy + 1;i++)
        delete newimg[i];
    delete newimg;
    newimg = NULL;
	// extract fea & normalize
	double***hog = new double**[sampleLen];
	for (int i = 0; i < sampleLen; i++)
	{
		hog[i] = new double*[ftrLen];
		for (int j = 0; j < ftrLen; j++)
		{
			hog[i][j] = new double[36];
			int x = ftrposx[j] + samplex[i] - imgminx-2;
			int y = ftrposy[j] + sampley[i] - imgminy-2;
			int w = ftrposw[j];
			int h = ftrposh[j];
            
			double norm = 0.00001;
			for (int k = 0; k < 9; k++)
			{
                hog[i][j][k] = iH[y + h/2-1][x + w/2-1][k];
                if(x>0&&y>0)
                    hog[i][j][k] = hog[i][j][k] - iH[y-1][x + w/2-1][k] - iH[y + h/2-1][x-1][k] + iH[y-1][x-1][k];
                else if(x>0)
                    hog[i][j][k] = hog[i][j][k] - iH[y + h/2-1][x-1][k];
                else if(y>0)
                    hog[i][j][k] = hog[i][j][k] - iH[y-1][x + w/2-1][k];
                norm += hog[i][j][k] * hog[i][j][k];
			}
            for (int k = 0; k < 9; k++)
			{
                hog[i][j][k+9] = iH[y + h/2-1][x + w-1][k] - iH[y+h/2-1][x + w/2-1][k];
                if(y>0)
                    hog[i][j][k+9] = hog[i][j][k+9]  - iH[y-1][x+w-1][k] + iH[y-1][x+w/2-1][k];
                norm += hog[i][j][k+9] * hog[i][j][k+9];
			}
             for (int k = 0; k < 9; k++)
			{
                hog[i][j][k+18] = iH[y + h-1][x + w/2-1][k] - iH[y+h/2-1][x + w/2-1][k];
                if(x>0)
                    hog[i][j][k+18] = hog[i][j][k+18] - iH[y + h-1][x-1][k] + iH[y+h/2-1][x-1][k];
                norm += hog[i][j][k+18] * hog[i][j][k+18];
			}
            for (int k = 0; k < 9; k++)
			{
                hog[i][j][k+27] = iH[y + h-1][x + w-1][k] - iH[y+h-1][x + w/2-1][k] - iH[y + h/2-1][x+w-1][k] + iH[y+h/2-1][x+w/2-1][k];
                norm += hog[i][j][k+27] * hog[i][j][k+27];
			}
			norm = sqrt(norm);
			for (int k = 0; k < 36; k++)
				hog[i][j][k] = hog[i][j][k]/norm;
        }
	}
        
	for (int i = 0; i < ftrLen; i++)
		for (int j = 0; j < sampleLen; j++)
			for (int k = 0; k < 36; k++)
				hog0[i + j*ftrLen + k*sampleLen*ftrLen] = hog[j][i][k];
    
    for (int i = 0; i < sampleLen; i++)
	{		
		for (int j = 0; j < ftrLen; j++)
            delete hog[i][j];
        delete hog[i];
    }
    delete hog;
    hog = NULL;
    for (int i = 0; i < imgmaxy - imgminy + 1; i++)
	{
		for (int j = 0; j < imgmaxx - imgminx + 1; j++)
			delete iH[i][j];
        delete iH[i];
    }
    delete iH;
    iH = NULL;
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	//int img0[], int w, int h, int samples0[], int sampleLen, int ftrpos0[], int ftrLen, double* hog0
	double * img0, *samples0, *ftrpos0;
	double *hog0;

	mwSize len_w, len_h, len_sampleLen, len_ftrLen;
	img0 = mxGetPr(prhs[0]);
	samples0 = mxGetPr(prhs[1]);
	ftrpos0 = mxGetPr(prhs[2]);
	
	len_h = mxGetM(prhs[0]);
	len_w = mxGetN(prhs[0]);
	len_sampleLen = mxGetM(prhs[1]);
	len_ftrLen = mxGetM(prhs[2]);


	plhs[0] = mxCreateDoubleMatrix(len_ftrLen*len_sampleLen, 36, mxREAL);

	hog0 = mxGetPr(plhs[0]);

    int* img = new int[len_h*len_w];
    for(int i=0;i<len_w;i++)
        for(int j=0;j<len_h;j++)
            img[i+j*len_w] = img0[i*len_h+j];
    int* samples = new int[len_sampleLen*2];
    for(int i=0;i<len_sampleLen*2;i++)
        samples[i] = samples0[i];
    int* ftrpos = new int[len_ftrLen*4];
    for(int i=0;i<len_ftrLen*4;i++)
        ftrpos[i] = ftrpos0[i];
	getFtrValHOG(img, len_w, len_h, samples, len_sampleLen, ftrpos, len_ftrLen, hog0);
    delete img;
    img = NULL;
    delete samples;
    samples = NULL;
    delete ftrpos;
    ftrpos = NULL;
	return;

}