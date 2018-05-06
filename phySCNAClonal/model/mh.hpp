#include<vector>
#include<cstring>
#include <string>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.hpp"
#include "cngenotype.hpp"
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

	int tReadNum;
	int nReadNum;

	bool baselineLabel;

	int copyNumber; //用于保存param时刻对应的copyNumber
	//free it before delete struct
	string genotype; //用于保存param时刻对应的genotype

	//此处不使用tp因为默认只使用一个tp
	double log_ll(double phi, cngenotype& cgn, int maxCopyNumber, double
			baseline){
		//pi 为基因型
		if(baselineLabel){
			ArrayXd cns(3);
			cns << 1, 2, 3;
			return log_likelihood_RD_BAF(phi, cgn, cns, baseline);
		}else if(get_loga(tReadNum, nReadNum)
				> baseline){
			ArrayXd cns(maxCopyNumber);
			for(int i = 1; i < maxCopyNumber; i ++){
				cns << i;
			}
			return log_likelihood_RD_BAF(phi, cgn, cns, baseline);
		}else{
			ArrayXd cns(3);
			cns << 0, 1, 2;
			return log_likelihood_RD_BAF(phi, cgn, cns, baseline);
		}
	}
	double log_likelihood_RD_BAF(double phi, cngenotype& cgn, ArrayXd& cns,
			double baseline){
		/************************************
		*  Here requires vector maximum operation  *
		************************************/
		int cnsLength = cns.size();
		int colsN = cns(cnsLength -1) + 1;
		int rowsN = cnsLength;

		ArrayXd barC = phi * cns + (1.0 - phi) * 2.0;
		ArrayXd lambdaPossion = (barC/2.0)*baseline*
			(nReadNum + 1);
		ArrayXd rds = log_poisson_pdf(lambdaPossion);

		/***************************************
		*  initialize with negative infinity  *
		***************************************/

		ArrayXXd ll(rowsN, colsN);
		ll.fill(-std::numeric_limits<double>::infinity());

		/********************************
		*  filter out b and d vector.  *
		********************************/

		//column vector
		ArrayXd bTj = get_b_T_j(a, b);
		ArrayXd dTj = get_d_T_j(a, b);

		for(int i=0; i<rowsN; i++){
			//muE should be row vector
			//mu_N constant
			//mu_G column vector
			//c_N constant
			//copyNumber constant
			//phi constant
			//ArrayXd muE = get_mu_E_joint(mu_N, mu_G, c_N, copyNumber, phi);
			ArrayXd muE = get_mu_E_joint(cgn.getBaf(cns(i)), copyNumber, phi);
			ll.block(i,0,1,cns(i)+1) = log_binomial_likelihood(bTj, dTj, muE);
		}

		/******************************
		*  get maximum args from ll  *
		******************************/

		MatrixXf::Index maxRow, maxCol;
		float maxLL = ll.maxCoeff(&maxRow, &maxCol);

		copyNumber = cns(maxRow);
		genotype = cgn.getGenotype(copyNumber, maxCol);

		return maxLL;
	}
};
