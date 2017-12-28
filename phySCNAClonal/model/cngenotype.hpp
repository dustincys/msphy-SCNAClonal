#ifndef CNGENOTYPE_H
#define CNGENOTYPE_H

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/SpecialFunctions>

using namespace Eigen;
using namespace std;

class cngenotype
{
public:
	cngenotype (int max_copy_number);

	string getGenotype(int copy_number, int index);
	vector<string> getGenotype(int copy_number);

	double getBaf(int copy_number, int index);
	ArrayXd getBaf(int copy_number);

private:
	//2: ["PP", "PM", "MM"]
	map< int, vector<string> > cn_genotype;
	//2: [0.2, 0.5, 0.2]
	//是否添加基因型先验概率表
	//map< int, vector<double> > cn_genotype_ll;
	//2: [0.2, 0.5, 0.2]
	map< int, ArrayXd> cn_genotype_baf;

	void init_map(int max_copy_number);
	string genotype(int p_num, int m_num);
};

string cngenotype::getGenotype(int copy_number, int index){
	return this->cn_genotype[copy_number][index];
}

vector<string> cngenotype::getGenotype(int copy_number){
	return this->cn_genotype[copy_number];
}

double cngenotype::getBaf(int copy_number, int index){
	return this->cn_genotype_baf[copy_number](index);
}

ArrayXd cngenotype::getBaf(int copy_number){
	return this->cn_genotype_baf[copy_number];
}

cngenotype::cngenotype(int max_copy_number){
	this->init_map(max_copy_number);
}

string cngenotype::genotype(int P_num, int M_num){
	string P_temp(P_num, 'P');
	string M_temp(M_num, 'M');
	return P_temp + M_temp;
}

void cngenotype::init_map(int max_copy_number){
	double EMPIRI_BAF =  0.485;
	double EMPIRI_AAF = 1 - EMPIRI_BAF;

	for(int cn=0; cn<max_copy_number+1; cn++){
		int rows_n = ceil((max_copy_number+1.0)/2.0);
		this->cn_genotype_baf[cn] = ArrayXd(rows_n);
		for(int M_num=0; M_num<(cn + 2)/2; M_num++){
			double mu_T = 0.0;
			string pi_T = "";

			int P_num = cn - M_num;
			if(0 == P_num && 0 == M_num){
				mu_T = EMPIRI_BAF/(EMPIRI_AAF + EMPIRI_BAF);
				pi_T = "NULL";
			}else if(P_num == M_num){
				mu_T = 0.5;
				pi_T = this->genotype(P_num, M_num);
			}else{
				mu_T = (M_num * 1.0) / cn;
				pi_T = this->genotype(P_num, M_num) +
					'/' + this->genotype(M_num, P_num);
			}
			this->cn_genotype[cn].push_back(pi_T);
			this->cn_genotype_baf[cn](M_num) = mu_T;
		}
	}
}
#endif
