#ifndef CNGENOTYPE_H
#define CNGENOTYPE_H 1

#include <string>
#include <iostream>
#include <vector>
#include <map>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/SpecialFunctions>

using namespace Eigen;
using namespace std;

class CNGenotype
{
public:
	CNGenotype (int maxCopyNumber);

	string getGenotype(int copyNumber, int index);
	vector<string> getGenotype(int copyNumber);

	double getBAF(int copyNumber, int index);
	ArrayXd getBAF(int copyNumber);

private:
	//2: ["PP", "PM", "MM"]
	map< int, vector<string> > cnGenotype;
	//2: [0.2, 0.5, 0.2]
	//是否添加基因型先验概率表
	//map< int, vector<double> > cn_genotype_ll;
	//2: [0.2, 0.5, 0.2]
	map< int, ArrayXd> cnGenotypeBAF;

	void init_map(int maxCopyNumber);
	string genotype(int pNum, int mNum);
};
#endif
