#ifndef MH_H
#define MH_H 1
#include<vector>
#include<cstring>
#include <string>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.hpp"
#include "cngenotype.hpp"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand);
double multi_param_post(struct node nodes[], struct datum data[], int old, struct config conf, CNGenotype& cngenotype);
double param_post(struct node nodes[], struct datum data[], int old, struct config conf, CNGenotype& cngenotype);
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
	int MAX_COPY_NUMBER; // maximum copy number
	float BASELINE; // baseline
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

class Stripe
{
public:

	string id;

	ArrayXd a;
	ArrayXd b;
	ArrayXXd bd;

	int tReadNum;
	int nReadNum;

	string tag; //时间戳标签

	int copyNumber; //用于保存param时刻对应的copyNumber
	string genotype; //用于保存param时刻对应的genotype



	Stripe (){};

	double log_ll(double phi, CNGenotype& cgn, int maxCopyNumber,
			double baseline){
		//pi 为基因型
		vector<int> cns;
		if(tag == "BASELINE"){
			for(int i=1: i<=3; i++){cns.push_back(i);}
		}else if(get_loga(this->tReadNum, this->nReadNum) > baseline){
			for(int i=2: i<=maxCopyNumber+1; i++){cns.push_back(i);}
		}else{
			for(int i=0: i<=2; i++){cns.push_back(i);}
		}

		ArrayXd lls(cns.size());
		ArrayXd gtIdxMaxs(cns.size());

		for(int i=0; i<=cns.size(); i++){
			lls[i] << this->getLLStripe(cns[i], phi, baseline, cgn,
					&gtIdxMaxs[i]);
		}

		int idxMax;
		double  ll = lls.maxCoeff(&idxMax);

		this->copyNumber = cns[idxMax];
		this->genotype = cgn.getGenotype(this->copyNumber,
				gtIdxMaxs[idxMax]);

		return ll;
	}

private:
	double getLLStripe(int copyNumber, double phi, double baseline,
			CNGenotype& cgn, int* gtIdxMax){
		double rdWeight = constants.RD_WEIGHT;
		double llStripe = 0;
		double llRD = this->getRD(copyNumber, phi, baseline);

		//此处在生成数据时使用numpy进行过滤即可，这里不进行操作
		//this->augBAF(copyNumber);

		ArrayXd bTj = get_b_T_j(this->a, this->b);
		ArrayXd dTj = get_d_T_j(this->a, this->b);

		double llBAF = 0;

		if(0 == bTj.size()){
			llBAF = 0;
			*gtIdxMax = -1;
		}else{
			llBAF = this->getBAF(phi, copyNumber, cgn, bTj, dTj,
					gtIdxMax);
		}

		double llRD = this->getRD(phi, copyNumber, baseline);

		return llRd * rdWeight + llBAF * (1 - rdWeight);
	}


	double getBAF(double phi, int copyNumber, CNGenotype& cgn,
			ArrayXd b, ArrayXd d, int *gtIdxMax){
		//muE row vector
		//size(b_T_j) x1
		//b_T_j column vector
		//d_T_j column vector
		//
		//return row vector
		//1 x size(muE)
		ArrayXd muE = get_mu_E_joint(cgn.getBAF(copyNumber),
				constants.MU_N, constants.COPY_NUMBER_NORMAL,
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

		float llBAF = llBAFs.maxCoeff(gtIdxMax);

		return llBAF;
	}

	double getRD(double phi, int copyNumber, double baseline){
		int cN = constants.COPY_NUMBER_NORMAL;
		double cMIN = constants.MINIMUM_POSITIVE;
		double barC = phi * copyNumber + (1.0 - phi) * cN;
		double lambdaPossion = (barC / cN) * baseline * (this->nReadNum
				+ 1.0);
		if(lambdaPossion <= 0){
			lambdaPossion = cMIN;
		}
		double llRD = log_poisson_pdf(this->tReadNum, lambdaPossion);
		return llRD;
	}
};
#endif
