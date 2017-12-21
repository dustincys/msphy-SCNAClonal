#include<vector>
#include<cstring>
#include <string>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.hpp"

using namespace std;

void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand,int tp);
double multi_param_post(struct node nodes[], struct datum data[], int old,struct config conf);
double param_post(struct node nodes[], struct datum data[], int old,struct config conf, int tp);
void update_params(struct node nodes[], struct config conf);
void get_pi(struct node nodes[], double pi[], struct config conf, int old, int tp);

void load_ssm_data(char fname[], struct datum data[], struct config conf);
void load_cnv_data(char fname[], struct datum data[], struct config conf);
void load_data_states(char fname[], struct datum data[], struct node nodes[], struct config conf);

void load_tree(char fname[], struct node nodes[], struct config conf);
void write_params(char fname[], struct node nodes[], struct config conf);

void mh_loop(struct node nodes[], struct datum data[], char* fname, struct config conf);

struct config{
	int MH_ITR;
	float MH_STD;

	int N_SSM_DATA; // no. of data points
	int N_CNV_DATA; // no. of data points

	int NNODES; // no. of nodes in the tree
	int TREE_HEIGHT;
	int NTPS; // no. of samples / time points
};

struct state{
	struct node* nd;
	int nr,nv;
};

struct node{
	int id;
	vector<double> param,pi;
	vector<double> param1,pi1; // dummy
	int ndata;
	vector<int> dids;
	int nchild;
	vector<int> cids; // children ids
	int ht;
};


struct datum{
	int id;

	//用于保存原始seg的索引，用于post process
	vector<int> segs_idx;

	vector<int> a;
	vector<int> b;

	int tumor_reads_num;
	int normal_reads_num;

	bool baseline_label;

	int copy_number; //用于保存param时刻对应的copy_number
	//free it before delete struct
	string genotype; //用于保存param时刻对应的genotype

	//vector<double> log_bin_norm_const;//log_bin_coeff(d,a);
	//struct datum* cnv; // for SSM datum, this is a pointer to its CNV datum
	//int cnv;// just an indicator for cnv or ssm datum
	// this is used to compute the binomial parameter
	//vector <struct state> states1, states2, states3, states4; // maternal and paternal state
	//double log_ll1111(vector<double> phi, int old){
	//double llh = 0.0;
	//for(int tp=0; tp<phi.size();tp++)
	//llh+=log_complete_ll(phi[tp],mu_r,mu_v,old,tp);
	//return llh;
	//}

	//此处不使用tp因为默认只使用一个tp
	double log_ll(double phi, cngenotype& cgn){
		//pi 为基因型
		ll, cn, pi = log_likelihood_RD_BAF(phi, cgn);

		copy_number = cn;
		genotype = pi;
		return ll
	}
	double log_likelihood_RD_BAF(double phi, cngenotype& cgn){
		/************************************
		*  Here requires metrix operation  *
		************************************/

        copy_numbers = None
        if seg.baseline_label == "True":
            copy_numbers = [2]
        elif get_loga(seg) > self._baseline:
            copy_numbers = range(2, self._max_copy_number + 1)
        else:
            copy_numbers = range(0, 2 + 1)

		return llh;

	}
};
