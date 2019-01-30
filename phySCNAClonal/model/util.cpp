#include <iostream>
#include <fstream>
#include <math.h>
#include "util.hpp"
#include "constants.hpp"

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
			+ 1.0);
}
double log_poisson_pdf(int tReadNum, double lambdaPossion){
	return tReadNum * log(lambdaPossion) - lambdaPossion -
		lgamma(tReadNum + 1.0);
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
	ArrayXXd ab(a.size(), 2);
	ab << a, b;
	return ab.rowwise().minCoeff();
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
        int cN = CONSTANTS::COPY_NUMBER_NORMAL;
        double muN = CONSTANTS::MU_N;
	ArrayXd muET = ((1.0 - phi) * cN * muN + phi * copyNumber * muG) /
		((1.0 - phi) * cN + phi * copyNumber);
	//muE should be row vector
	return muET.transpose();
}

void getLLStripe(int segIdxSize, int tag, int copyNumber, double phi, double
		baseline, CNGenotype cgn, int& gtIdxMax, ArrayXd a, ArrayXd b,
		int tReadNum, int nReadNum, double& ll){
	double rdWeight = CONSTANTS::RD_WEIGHT;
	double llStripe = 0;
	double llRD = getRD(copyNumber, phi, baseline, tReadNum, nReadNum);

	//此处在生成数据时使用numpy进行过滤即可，这里不进行操作
	//augBAF(copyNumber);

	ArrayXd bTj = get_b_T_j(a, b);
	ArrayXd dTj = get_d_T_j(a, b);

	double llBAF = 0;

	if(0 == bTj.size()){
		llBAF = 0;
		gtIdxMax = -1;
	}else{
		llBAF = getBAF(phi, copyNumber, cgn, bTj, dTj,
				gtIdxMax);
	}

	//cout << "llRD" << llRD << endl;
	//cout << "llBAF" << llBAF << endl;
	//cout << "rdWeight" << rdWeight << endl;
	ll = llRD * rdWeight * segIdxSize+ llBAF * (1 - rdWeight);
	return;
}


double getBAF(double phi, int copyNumber, CNGenotype cgn,
		ArrayXd b, ArrayXd d, int& gtIdxMax){
	//muE row vector
	//size(b_T_j) x1
	//b_T_j column vector
	//d_T_j column vector
	//
	//return row vector
	//1 x size(muE)
	ArrayXd muE = get_mu_E_joint(cgn.getBAF(copyNumber),
			CONSTANTS::MU_N, CONSTANTS::COPY_NUMBER_NORMAL,
			copyNumber, phi);
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
		(dArray - bArray + 1).lgamma() + bArray * muArray.log()
		+ (dArray - bArray) * (1 - muArray).log();

	/*--  Here returns CN ll vector and the best genotype vector
	 * --*/
	ArrayXd llBAFs = ll.matrix().colwise().sum();
	float llBAF = llBAFs.maxCoeff(&gtIdxMax);

	return llBAF;
}

double getRD(int copyNumber, double phi, double baseline, int tReadNum, int
		nReadNum){
	int cN = CONSTANTS::COPY_NUMBER_NORMAL;
	double cMIN = CONSTANTS::MINIMUM_POSITIVE;
	double barC = phi * copyNumber + (1.0 - phi) * cN;
	double lambdaPossion = (barC / cN) * exp(baseline) * (nReadNum +
			1.0);
	if(lambdaPossion <= 0){
		lambdaPossion = cMIN;
	}
	double llRD = log_poisson_pdf(tReadNum, lambdaPossion);
	return llRD;
}

