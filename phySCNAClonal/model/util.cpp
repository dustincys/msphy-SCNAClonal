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

void dirichlet_sample(int size, double alpha[], double *x,gsl_rng *r){
	gsl_ran_dirichlet(r,size,alpha,x);
	for(int i=0;i<size;i++)
		x[i]=x[i]+0.0001;
	//此处是非零处理
	double sum=0.0;
	for(int i=0;i<size;i++)
		sum+=x[i];
	for(int i=0;i<size;i++)
		x[i]/=sum;
}

double logsumexp(double x[], int nx){
	double maxes=x[0], sum=0.0;
	for (int i=1;i<nx;i++)
		if (x[i]>maxes)
			maxes=x[i];
	for (int i=0;i<nx;i++)
		sum+=exp(x[i]-maxes);
	return log(sum) + maxes;
}

vector<int> get_cn_vec(int max_copy_number){
	vector<int> cn_vec;
	for(int i=0; i<max_copy_number+1; i++){
		cn_vec.push_back(i)
	}
	return cn_vec;
}

double get_loga(int tumor_reads_num, int normal_reads_num){
	if(0 == tumor_reads_num){
		tumor_reads_num = 1;
	}
	if(0 == normal_reads_num){
		normal_reads_num = 1;
	}

	return log(tumor_reads_num) - log(normal_reads_num);
}

ArrayXd log_poisson_pdf(ArrayXd lambda_possion){
	return tumor_reads_num *lambda_possion.log() - lambda_possion
		- lgamma(tumor_reads_num + 1);
}

ArrayXd get_mu_E_joint(mu_N, mu_G, c_N, c_H, phi){
	return ((1.0 - phi) * c_N * mu_N + phi * c_H * mu_G) / (
			(1.0 - phi) * c_N + phi * c_H);
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

//mu_N constant
//mu_G column vector
//c_N constant
//copy_number constant
//phi constant
//ArrayXd mu_E = get_mu_E_joint(mu_N, mu_G, c_N, copy_number, phi);
ArrayXd get_mu_E_joint(ArrayXd mu_G, int copy_number, double phi){
        int c_N = constants.COPY_NUMBER_NORMAL;
        double mu_N = constants.MU_N;
	ArrayXd mu_E_T = ((1.0 - phi)*c_N*mu_N + phi*copy_number*mu_G)/((1.0 -
				phi)*c_N + phi*copy_number);
	//mu_E should be row vector
	return mu_E_T.transpose();
}


ArrayXd log_binomial_likelihood(ArrayXd b, ArrayXd d, ArrayXd mu_E){
	//mu_E row vector
	//b_T_j column vector
	//d_T_j column vector
	//
	//return row vector

	//1 x size(mu_E)
	MatrixXd v1 = (ArrayXXd::Zero(1,3) + 1).matrix();
	//n x 1
	ArrayXXd b_array = (b.matrix() * v1).array();
	ArrayXXd d_array = (d.matrix() * v1).array();
	ArrayXXd mu_array = (v1.transpose() * mu_E).array();

	ArrayXXd ll = (d_array + 1).lgamma() - (b_array + 1).lgamma() -
		(d_array - b_array).lgamma() + b_array * mu_array.log() +
		(d_array - b_array) * (1 - mu_array);

	return ll.colwise().sum();
}
