#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <map>
#include <algorithm>


#include "mh.hpp"
#include "util.hpp"

using namespace std;

//  g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`
// ./mh.o 5000 100 11 1 8 4 ssm_data.txt cnv_data.txt c_tree_ssm_data_1.txt c_data_states_ssm_data_1.txt c_params.txt 5
//https://www.gnu.org/software/gsl/manual/html_node/Shared-Libraries.html


// done for multi-sample
int main(int argc, char* argv[]){
	// parse command line args
	struct config conf;
	conf.MH_ITR=atoi(argv[1]);//5000 iteration
	conf.MH_STD=atof(argv[2]);//100 STD
	conf.N_STRIPE_DATA=atoi(argv[3]);//12; // no. of stripe data points
	conf.NNODES=atoi(argv[4]); //17, no. of nodes in the tree
	conf.TREE_HEIGHT=atoi(argv[5]);//6

	char* FNAME_STRIPE_DATA = argv[6];
	char* FNAME_C_TREE = argv[7];
	char* FNAME_C_PARAMS = argv[8];
	char* FNAME_C_MH_AR = argv[9];

	struct datum *data = new datum[conf.N_STRIPE_DATA];
	load_stripe_data(FNAME_STRIPE_DATA, data, conf);

	struct node *nodes = new node[conf.NNODES];
	load_tree(FNAME_C_TREE,nodes,conf);

	//start MH loop
	mh_loop(nodes,data,FNAME_C_MH_AR,conf);

	// write updated params to disk
	write_params(FNAME_C_PARAMS,nodes,conf);

	return 0;
}


// done for multi-sample
void mh_loop(struct node nodes[],struct datum data[],char* fname, struct config conf){
	gsl_rng *rand = gsl_rng_alloc(gsl_rng_mt19937);
	double ratio=0.0;
	for (int itr=0;itr<conf.MH_ITR;itr++){

		sample_cons_params(nodes,conf,rand);

		double a = multi_param_post(nodes,data,0,conf)-multi_param_post(nodes,data,1,conf);

		//cout<<multi_param_post(nodes,data,0,conf)-multi_param_post(nodes,data,1,conf)<<'\n';
		// loop over samples, apply dirichlet correction terms, update a
		double theta[conf.NNODES];// dirichlet params
		double pi[conf.NNODES],pi_new[conf.NNODES];

		get_pi(nodes,pi_new,conf,0);
		get_pi(nodes,pi,conf,1);

		// apply the dirichlet correction terms
		for(int i=0;i<conf.NNODES;i++)
			theta[i]=conf.MH_STD*pi_new[i];
		a += gsl_ran_dirichlet_lnpdf(conf.NNODES,theta,pi);

		for(int i=0;i<conf.NNODES;i++)
			theta[i]=conf.MH_STD*pi[i];
		a -= gsl_ran_dirichlet_lnpdf(conf.NNODES,theta,pi_new);

		double r = gsl_rng_uniform_pos(rand);
		//cout<<log(r)<<'\t'<<a<<'\n';
		if (log(r)<a){
			ratio+=1;
			update_params(nodes,conf);
		}
	}
	gsl_rng_free(rand);

	ofstream dfile;
	dfile.open(fname);
	dfile<<ratio/conf.MH_ITR;
	dfile.close();
}


// done for multi-sample
void sample_cons_params(struct node nodes[], struct config conf, gsl_rng *rand){

	map <int, int> node_id_map;
	for(int i=0;i<conf.NNODES;i++)
		node_id_map[nodes[i].id]=i;

	int NNODES=conf.NNODES;
	double pi[NNODES];
	for(int i=0;i<NNODES;i++)
		pi[i]=nodes[i].pi;

	// randomly sample from a dirichlet
	double pi_new[NNODES],alpha[NNODES];
	for(int i=0;i<NNODES;i++)
		alpha[i]=conf.MH_STD*pi[i]+1;
	dirichlet_sample(NNODES,alpha,pi_new,rand);

	// update the nodes pi1 (new pi)
	for(int i=0;i<NNODES;i++)
		nodes[i].pi1=pi_new[i];

	// update the nodes param1 (new param)
	for(int i=0;i<NNODES;i++){
		double param = nodes[i].pi1;
		for(int c=0;c<nodes[i].nchild;c++){
			param+=nodes[node_id_map[nodes[i].cids.at(c)]].param1;
		}
		nodes[i].param1=param;
	}
}


// done for multi-sample
// todo: double check log_ll
double multi_param_post(struct node nodes[], struct datum data[], int old, struct config conf){
	return param_post(nodes,data,old,conf);
}
double param_post(struct node nodes[], struct datum data[], int old,struct config conf){
	double llh = 0.0;
	for(int i=0;i<conf.NNODES;i++){
		double p=0;
		if(old==0)
			p=nodes[i].param1;
		else
			p=nodes[i].param;
		for(int j=0;j<nodes[i].ndata;j++){
			//此处调用data中的似然 data id
			llh+=data[nodes[i].dids.at(j)].log_ll(p, old);
		}
	}
	return llh;
}


// done for multi-sample
void update_params(struct node nodes[],struct config conf){
	for(int i=0;i<conf.NNODES;i++){
		nodes[i].param=nodes[i].param1;
		nodes[i].pi=nodes[i].pi1;
	}
}

// done for multi-sample
void get_pi(struct node nodes[], double pi[], struct config conf, int old){
	for(int i=0;i<conf.NNODES;i++){
		if (old==0)
			pi[i]=nodes[i].pi1;
		else
			pi[i]=nodes[i].pi;
	}
}


// done for multi-sample
void write_params(char fname[], struct node *nodes, struct config conf){
	ofstream dfile;
	dfile.open(fname);
	for(int i=0;i<conf.NNODES;i++){
		dfile<<nodes[i].id<<'\t';
		dfile<<nodes[i].param<<',';
		dfile<<'\t';
		dfile<<nodes[i].pi<<',';
		dfile<<'\n';
	}
	dfile.close();
}


// done for multi-sample
void load_stripe_data(char fname[],struct datum *data, struct config conf){
	string line,token,token1,token2;
	ifstream dfile (fname);
	int ctr=0,id=-1,ab=1;
	while (getline (dfile,line,'\n')){
		if (id==-1){id+=1;continue;}
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				data->id=id;
			}
			else if(ctr==2){
				istringstream iss(token);
				ab = 1;
				while(getline(iss,token1,'|')){
					istringstream iss1(token1);
					size_t n = (std::count(t.begin(), s.end(), ',') + 2)/2;
					if(ab==1){
						data->a = ArrayXd(n);
						data->a << atoi(token1.c_str());
					}else{
						data->b = ArrayXd(n);
						data->b << atoi(token1.c_str());
					}
					ab++;
				}
			}
			else if(ctr==3){
				istringstream iss(token);
				data->tumor_reads_num=atof(token.c_str());
			}
			else if(ctr==4){
				istringstream iss(token);
				data->normal_reads_num=atof(token.c_str());
			}
			else if(ctr==5){
				istringstream iss(token);
				data->baseline_label=atob(token.c_str());
			}
			ctr+=1;
		}
		data++;
		id+=1;
	}
	dfile.close();
}

void load_tree(char fname[], struct node *nodes, struct config conf){
	string line,token,token1;
	ifstream dfile (fname);
	int ctr=0;
	while (getline (dfile,line,'\n')){
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				nodes->id=atoi(token.c_str());
			}
			else if(ctr==1){
				nodes->param.push_back(atof(token.c_str()));
				nodes->param1.push_back(0.0);
			}
			else if(ctr==2){
				nodes->pi.push_back(atof(token.c_str()));
				nodes->pi1.push_back(0.0);
			}
			else if(ctr==5){
				nodes->ndata=atoi(token.c_str());
			}
			else if(ctr==6){
				istringstream iss(token);
				for(int i=0;i<nodes->ndata;i++){
					getline(iss,token1,',');
					nodes->dids.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==3){
				nodes->nchild=atoi(token.c_str());
			}
			else if(ctr==4){
				istringstream iss(token);
				for(int i=0;i<nodes->nchild;i++){
					getline(iss,token1,',');
					nodes->cids.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==7){
				nodes->ht=atoi(token.c_str());
			}
			ctr+=1;
		}
		nodes++;
	}
	dfile.close();
}
