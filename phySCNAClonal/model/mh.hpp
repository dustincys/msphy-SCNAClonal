#ifndef MH_H
#define MH_H 1
#include<vector>
#include<cstring>
#include <string>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.hpp"
//#include "cngenotype.hpp"
#include "constants.hpp"
#include "Eigen/Dense"
#include <thread>         // std::thread, std::this_thread::sleep_for

using namespace std;
using namespace Eigen;


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

class SCNA
{
public:

	string id;

	ArrayXd a;
	ArrayXd b;

	int tReadNum;
	int nReadNum;

	string tag; //时间戳标签

	int copyNumber; //用于保存param时刻对应的copyNumber
	string genotype; //用于保存param时刻对应的genotype



	SCNA (){};

	double log_ll(double phi, CNGenotype& cgn, int maxCopyNumber,
			double baseline){
		//pi 为基因型
		vector<int> cns;
		if(tag == "BASELINE"){
			for(int i=1; i<=3; i++){cns.push_back(i);}
		}else if(get_loga(this->tReadNum, this->nReadNum) > baseline){
			for(int i=2; i<=maxCopyNumber; i++){cns.push_back(i);}
		}else{
			for(int i=0; i<=2; i++){cns.push_back(i);}
		}

		double lls[cns.size()];
		int gtIdxMaxs[cns.size()];

		thread threads[cns.size()];                         // default-constructed threads
		for(int i=0; i<cns.size(); i++){
			threads[i] = thread(getLLStripe, cns[i], phi, baseline,
					cgn, ref(gtIdxMaxs[i]), this->a, this->b,
					this->tReadNum, this->nReadNum, ref(lls[i]));
		}
		for(int i=0; i<cns.size(); i++){
			threads[i].join();
		}

		double ll = *max_element(lls, lls+cns.size());
		int idxMax = distance(lls, max_element(lls, lls+cns.size()));
		//double  ll = lls.maxCoeff(&idxMax);

		this->copyNumber = cns[idxMax];

		if (gtIdxMaxs[idxMax] == -1) {
			this->genotype = "UNCERTAIN";
		}else{
			this->genotype = cgn.getGenotype(this->copyNumber,
					gtIdxMaxs[idxMax]);
		}

		return ll;
	}

};




void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand);
double multi_param_post(struct node nodes[], SCNA data[], int old, struct config conf, CNGenotype& cngenotype);
double param_post(struct node nodes[], SCNA data[], int old, struct config conf, CNGenotype& cngenotype);
void update_params(struct node nodes[], struct config conf);
void get_pi(struct node nodes[], double pi[], struct config conf, int old);
void load_stripe_data(char fname[],SCNA *data, struct config conf);
void load_tree(char fname[], struct node nodes[], struct config conf);
void write_params(char fname[], struct node nodes[], struct config conf);
void mh_loop(struct node nodes[], SCNA data[], char* fname, struct config conf);
#endif
