#include <iostream>
#include<math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#ifndef UTIL_H
#define UTIL_H 1

struct config_constant{
	int COPY_NUMBER_NORMAL = 2;
        double MU_N = 0.5;
}constants;

// stat/math functions
double log_factorial(int n);
double log_bin_coeff(int n, int k);
double log_binomial_likelihood(int x, int n, double mu);
double log_beta(double a, double b);
void dirichlet_sample(int size, double alpha[], double *x,gsl_rng *r);
double logsumexp(double x[], int nx);
double get_loga(int tumor_reads_num, int normal_reads_num);

vector<int> get_cn_vec(int max_copy_number);

ArrayXd get_b_T_j(ArrayXd a, ArrayXd b);
ArrayXd get_d_T_j(ArrayXd a, ArrayXd b);
ArrayXd get_mu_E_joint(ArrayXd mu_G, int copy_number, double phi);
ArrayXd log_binomial_likelihood(ArrayXd b, ArrayXd d, ArrayXd mu_E);
ArrayXd log_poisson_pdf(ArrayXd lambda_possion);

#endif
