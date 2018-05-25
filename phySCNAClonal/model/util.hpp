#ifndef UTIL_H
#define UTIL_H 1
#include <iostream>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include "Eigen/Dense"
#include "constants.hpp"


using namespace Eigen;
using namespace std;


// stat/math functions
double log_factorial(int n);
double log_bin_coeff(int n, int k);
double log_binomial_likelihood(int x, int n, double mu);
double log_beta(double a, double b);
void dirichlet_sample(int size, double alpha[], double *x, gsl_rng *r);
double logsumexp(double x[], int nx);
double get_loga(int tReadNum, int nReadNum);
vector<int> get_cn_vec(int maxCopyNumber);
ArrayXd get_b_T_j(ArrayXd a, ArrayXd b);
ArrayXd get_d_T_j(ArrayXd a, ArrayXd b);
ArrayXd get_mu_E_joint(ArrayXd muG, double muN, int cN, int cH, double phi);
ArrayXd log_binomial_likelihood(ArrayXd b, ArrayXd d, ArrayXd muE);
ArrayXd log_poisson_pdf(int tReadNum, ArrayXd lambdaPossion);
double log_poisson_pdf(int tReadNum, double lambdaPossion);
#endif
