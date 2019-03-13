#include <cmath>
#include <ctime>
#include "mex.h"
using namespace std;

class BPAClassifier
{
public:
	void train(double*beta, double**x, int*y, int N, int dim, double* w, double C);
	void predict(double**x, int N, double* pre);
	void getClassifier(double*newbeta);
    ~BPAClassifier();
private:
	void init(double*beta, double**x, int*y, int N, int dim, double* w, double C);
	bool satisfy_KKT(int id);
	double linearkernel(double* x1, double* x2);
	double get_error(int id);
	double getTrainLoss();
	int update(int x1);
	double*beta;
	double **x;
	int *y;
	int N;
	int dim;
	double* w;
	double C;
	double* alpha;
	double* newbeta;
	double** K;
	double* P;
 	const double TOL = 1e-5;
	const double EPS = 1e-5;
};
BPAClassifier::~BPAClassifier()
{
    delete beta;
    beta = NULL;
    for(int i=0;i<N;i++)
        delete x[i];
    delete x;
    x= NULL;
	delete y;
    y = NULL;
	
	delete w;
    w = NULL;
	delete alpha;
    alpha = NULL;
	delete newbeta;
    newbeta = NULL;
    for(int i=0;i<N;i++)
        delete K[i];
    delete K;
    K = NULL;
	delete P;
    P = NULL;
}
void BPAClassifier::getClassifier(double*newbeta)
{
	for (int i = 0; i < dim; i++)
		newbeta[i] = this->newbeta[i];
}
void BPAClassifier::predict(double**x, int N, double* pre)
{
	for (int i = 0; i < N; i++)
	{
		pre[i] = linearkernel(newbeta, x[i]);
	}
}

int BPAClassifier::update(int x1)
{
	double new_alpha = 1 - y[x1] * P[x1];
	for (int i = 0; i < N; i++)
		if (i != x1)
			new_alpha += -K[i][x1] * y[i] * y[x1] * alpha[i];
	new_alpha = new_alpha /  K[x1][x1];
	if (new_alpha > w[x1] * C) new_alpha = w[x1] * C;
	else if (new_alpha < 0) new_alpha = 0;
	//变化太小的话当做是没有发生变化
	if (abs(new_alpha - alpha[x1]) < EPS)
		return 0;
	alpha[x1] = new_alpha;
	return 1;
}
double BPAClassifier::linearkernel(double* x1, double* x2) {
	double result(0);
	for (int i = 0; i < dim; ++i) {
		result += (x1[i] * x2[i]);
	}
	return result;
}

double BPAClassifier::get_error(int id) {
	double rst = P[id];
	for (int i = 0; i < N; ++i) {
		if (alpha[i] != 0) {
			rst += (alpha[i] * y[i] * K[i][id]);
		}
	}
	return rst - y[id];
}

bool BPAClassifier::satisfy_KKT(int id) {
	double r = y[id] * get_error(id);
	if ((alpha[id]<C*w[id] && r < -TOL) || (alpha[id]>0 && r>TOL)) return false;
	return true;
}

void BPAClassifier::init(double*beta, double**x, int*y, int N, int dim, double* w, double C)
{
	this->dim = dim;
	this->N = N;
	this->w = new double[N];
	this->beta = new double[dim];
	this->x = new double*[N];
	this->y = new int[N];
	for (int i = 0; i < N; i++)
	{
		this->x[i] = new double[dim];
		for (int j = 0; j < dim; j++)
			this->x[i][j] = x[i][j];
		this->w[i] = w[i];
		this->y[i] = y[i];
	}
	for (int j = 0; j < dim; j++)
		this->beta[j] = beta[j];
	this->C = C;
	//计算线性核以及原始权重预测
	K = new double*[N];
	P = new double[N];
	for (int i = 0; i < N; i++)
	{
		K[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			if (i > j)
				K[i][j] = K[j][i];
			else
				K[i][j] = linearkernel(x[i], x[j]);
		}
		P[i] = linearkernel(beta, x[i]);
	}

	alpha = new double[N];
	for (int i = 0; i < N; i++)
		alpha[i] = 0;
}

void BPAClassifier::train(double*beta, double**x, int*y, int N, int dim, double* w, double C)
{

	init(beta, x, y, N, dim, w, C);

	newbeta = new double[dim];
	bool loop_on_bounds = false;
	int num_changed;		//没有直接用num_violate_KKT，因为有可能虽然违反了KKT，但是违反程度很小，a1 a2的变化非常细微。
	while (true) {
		num_changed = 0;
		for (int x1 = 0; x1 < N; ++x1) {
			if ((alpha[x1] == 0 || alpha[x1] == C*w[x1]) && loop_on_bounds) continue;	//只在边界上面寻找a
			if (!satisfy_KKT(x1)) 
				num_changed += update(x1);
		}
		if ((!loop_on_bounds) && num_changed == 0) break;	//在整个数据集上都找不到违反KKT条件的a，结束算法
		if (num_changed) loop_on_bounds = true;			//在整个数据集上遍历一次后，立刻回到边界点上遍历，直到边界点都符合KKT
		else			loop_on_bounds = false;
	}

	for (int i = 0; i < dim; i++)
	{
		newbeta[i] = beta[i];
		for (int j = 0; j < N; j++)
			newbeta[i] += alpha[j] * y[j] * x[j][i];
	}
}

double BPAClassifier::getTrainLoss()
{
	for (int i = 0; i < dim; i++)
	{
		newbeta[i] = beta[i];
		for (int j = 0; j < N; j++)
			newbeta[i] += alpha[j] * y[j] * x[j][i];
	}
	double* pre = new double[N];
	predict(x, N, pre);
	double loss = 0;
	for (int i = 0; i < N; i++)
		if (pre[i] * y[i] < 1)
			loss += (1 - pre[i] * y[i])*C*w[i];
	return loss;
}

static void train(double*beta, double**x, int*y, int N, int dim, double* w, double* cls)
{
    BPAClassifier model;
    model.train(beta, x, y, N, dim, w, 1);
    model.getClassifier(cls);   
}          
void mexFunction( int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] )

{ 
    if (nrhs!=4)
        mexErrMsgTxt("The number of parameters must be 4!\n");

    double*beta = mxGetPr(prhs[0]);
    double*tx = mxGetPr(prhs[1]);
    double*ty = mxGetPr(prhs[2]);
    double* w = mxGetPr(prhs[3]);
    

	size_t N, dim;

	N = mxGetM(prhs[1]);
	dim = mxGetN(prhs[1]);
    
    double**x = new double*[N];
    int*y = new int[N];
    for(int i=0;i<N;i++)
    {
        x[i] = new double[dim];
        for(int j=0;j<dim;j++)
            x[i][j] = tx[j*N+i];
        y[i] = ty[i];
    }
    
    plhs[0] = mxCreateDoubleMatrix(dim,1, mxREAL);
    double* wb = mxGetPr(plhs[0]);
	train(beta, x, y, N, dim, w, wb);
	
    for(int i=0;i<N;i++)
        delete x[i];
    delete x;
    x = NULL;
    delete y;
    y = NULL;
	return;

}