#include <iostream>
#include <fstream>
#include "util.hpp"

using namespace std;

// stat/math functions
double log_factorial(int n){
	return lgamma(n + 1);
}

double log_bin_coeff(int n, int k){
	return log_factorial(n) - log_factorial(k) - log_factorial(n - k);
}

double log_binomial_likelihood(int x, int n, double mu){
	return  x * log(mu) + (n - x) * log(1 - mu);
}

double log_beta(double a, double b){
	return lgamma(a) + lgamma(b) - lgamma(a + b);
}

void dirichlet_sample(int size, double alpha[], double *x, gsl_rng *r){
	gsl_ran_dirichlet(r, size, alpha, x);
	for(int i = 0; i < size; i++)
		x[i] = x[i] + 0.0001;
	//此处是非零处理
	double sum = 0.0;
	for(int i = 0; i < size; i++)
		sum += x[i];
	for(int i = 0; i < size; i++)
		x[i] /= sum;
}

double logsumexp(double x[], int nx){
	double maxes = x[0], sum = 0.0;
	for (int i = 1; i < nx; i++)
		if (x[i] > maxes)
			maxes = x[i];
	for (int i = 0; i < nx; i++)
		sum += exp(x[i] - maxes);
	return log(sum) + maxes;
}

vector<int> get_cn_vec(int maxCopyNumber){
	vector<int> cnVec;
	for(int i=0; i<maxCopyNumber+1; i++){
		cnVec.push_back(i);
	}
	return cnVec;
}

double get_loga(int tReadNum, int nReadNum){
	if(0 == tReadNum){
		tReadNum = 1;
	}
	if(0 == nReadNum){
		nReadNum = 1;
	}

	return log(tReadNum) - log(nReadNum);
}

ArrayXd log_poisson_pdf(int tReadNum, ArrayXd lambdaPossion){
	return tReadNum *lambdaPossion.log() - lambdaPossion - lgamma(tReadNum
			+ 1);
}

ArrayXd get_mu_E_joint(ArrayXd muG, double muN, int cN, int cH, double phi){
	return ((1.0 - phi) * cN * muN + phi * cH * muG) / ((1.0 - phi) * cN +
			phi * cH);
}

ArrayXd get_b_T_j(ArrayXd a, ArrayXd b){
	/******************************
	*  here, a, b column vector.  *
	******************************/
	assert(a.size() == b.size());
	ArrayXd ab(a.size(), 2);
	ab << a, b;
	return ab.rowwise().maxCoeff();
}

ArrayXd get_d_T_j(ArrayXd a, ArrayXd b){
	/******************************
	*  here, a, b column vector.  *
	******************************/
	return a + b;
}

//muN constant
//muG column vector
//cN constant
//copyNumber constant
//phi constant
//ArrayXd muE = get_mu_E_joint(muN, muG, cN, copyNumber, phi);
ArrayXd get_mu_E_joint(ArrayXd muG, int copyNumber, double phi){
        int cN = constants.COPY_NUMBER_NORMAL;
        double muN = constants.MU_N;
	ArrayXd muET = ((1.0 - phi) * cN * muN + phi * copyNumber * muG) /
		((1.0 - phi) * cN + phi * copyNumber);
	//muE should be row vector
	return muET.transpose();
}


ArrayXd log_binomial_likelihood(ArrayXd b, ArrayXd d, ArrayXd muE){
	//muE row vector
	//size(b_T_j) x1
	//b_T_j column vector
	//d_T_j column vector
	//
	//return row vector

	//1 x size(muE)
	MatrixXd v1 = (ArrayXXd::Zero(1, muE.size()) + 1).matrix();
	MatrixXd v2 = (ArrayXXd::Zero(b.size(), 1) + 1).matrix();
	////n x 1
	ArrayXXd bArray = (b.matrix() * v1).array();
	ArrayXXd dArray = (d.matrix() * v1).array();

	////此处需要注意维度匹配
	////muE Nx1
	////v1 1xN

	ArrayXXd muArray = (v2 * muE.transpose().matrix()).array();

	ArrayXXd ll = (dArray + 1).lgamma() - (bArray + 1).lgamma() -
		(dArray - bArray + 1).lgamma() + bArray * muArray.log() +
		(dArray - bArray) * (1 - muArray).log();


	return ll.matrix().colwise().sum();
}
