#include<vector>
#include<cstring>
#include <string>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.hpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand);
double multi_param_post(struct node nodes[], struct datum data[], int old,struct config conf);
double param_post(struct node nodes[], struct datum data[], int old,struct config conf);
void update_params(struct node nodes[], struct config conf);
void get_pi(struct node nodes[], double pi[], struct config conf, int old);

void load_stripe_data(char fname[],struct datum *data, struct config conf);

void load_tree(char fname[], struct node nodes[], struct config conf);
void write_params(char fname[], struct node nodes[], struct config conf);

void mh_loop(struct node nodes[], struct datum data[], char* fname, struct config conf);

struct config{
	int MH_ITR;
	float MH_STD;
	int N_STRIPE_DATA; // no. of data points
	int NNODES; // no. of nodes in the tree
	int TREE_HEIGHT;
};

struct node{
	int id;
	double param,pi;
	double param1,pi1; // dummy
	int ndata;
	vector<int> dids;
	int nchild;
	vector<int> cids; // children ids
	int ht;
};


struct datum{
	int id;

	//用于保存原始seg的索引，用于post process
	//vector<int> segs_idx;

	ArrayXd a;
	ArrayXd b;

	int tumor_reads_num;
	int normal_reads_num;

	bool baseline_label;

	int copy_number; //用于保存param时刻对应的copy_number
	//free it before delete struct
	string genotype; //用于保存param时刻对应的genotype

	//此处不使用tp因为默认只使用一个tp
	double log_ll(double phi, cngenotype& cgn,
			int max_copy_number, double baseline){
		//pi 为基因型
		if(baseline_label){
			ArrayXd cns(3);
			cns << 1, 2, 3;
			log_likelihood_RD_BAF(phi, cgn, cns);
		}else if(get_loga(tumor_reads_num, normal_reads_num)
				> baseline){
			ArrayXd cns(max_copy_number);
			for(int i = 1; i < max_copy_number; i ++){
				cns << i;
			}
			log_likelihood_RD_BAF(phi, cgn, cns);
		}else{
			ArrayXd cns(3);
			cns << 0, 1, 2;
			log_likelihood_RD_BAF(phi, cgn, cns);
		}
	}
	double log_likelihood_RD_BAF(double phi, cngenotype& cgn, ArrayXd& cns){
		/************************************
		*  Here requires vector maximum operation  *
		************************************/
		int cns_length = cns.size();
		int cols_n = cns(cns_length -1) + 1;
		int rows_n = cns_length;

		ArrayXd bar_c = phi * cns + (1.0 - phi) * 2.0;
		ArrayXd lambda_possion = (bar_c/2.0)*baseline*
			(normal_reads_num + 1);
		ArrayXd rds = log_poisson_pdf(lambda_possion);

		/***************************************
		*  initialize with negative infinity  *
		***************************************/

		ArrayXXd ll(rows_n, cols_n);
		ll.fill(-std::numeric_limits<double>::infinity());

		/********************************
		*  filter out b and d vector.  *
		********************************/

		//column vector
		ArrayXd b_T_j = get_b_T_j(a, b);
		ArrayXd d_T_j = get_d_T_j(a, b);

		for(int i=0; i<rows_n; i++){
			//mu_E should be row vector
			//mu_N constant
			//mu_G column vector
			//c_N constant
			//copy_number constant
			//phi constant
			//ArrayXd mu_E = get_mu_E_joint(mu_N, mu_G, c_N, copy_number, phi);
			ArrayXd mu_E = get_mu_E_joint(cgn.getBaf(cns(i)), copy_number, phi);
			ll.block(i,0,1,cns(i)+1) = log_binomial_likelihood(b_T_j, d_T_j, mu_E);
		}

		/******************************
		*  get maximum args from ll  *
		******************************/

		MatrixXf::Index maxRow, maxCol;
		float max_ll = ll.maxCoeff(&maxRow, &maxCol);

		copy_number = cns(maxRow);
		genotype = cgn.getGenotype(copy_number, maxCol);

		return max_ll;
	}
};
