#include "cngenotype.hpp"

string CNGenotype::getGenotype(int copyNumber, int index){
	return this->cnGenotype[copyNumber][index];
}

vector<string> CNGenotype::getGenotype(int copyNumber){
	return this->cnGenotype[copyNumber];
}

double CNGenotype::getBAF(int copyNumber, int index){
	return this->cnGenotypeBAF[copyNumber](index);
}

ArrayXd CNGenotype::getBAF(int copyNumber){
	return this->cnGenotypeBAF[copyNumber];
}

CNGenotype::CNGenotype(int maxCopyNumber){
	this->init_map(maxCopyNumber);
}

string CNGenotype::genotype(int pNum, int mNum){
	string pTemp(pNum, 'P');
	string mTemp(mNum, 'M');
	return pTemp + mTemp;
}

void CNGenotype::init_map(int maxCopyNumber){
	double EMPIRI_BAF =  0.485;
	double EMPIRI_AAF = 1 - EMPIRI_BAF;

	for(int cn=0; cn<maxCopyNumber+1; cn++){
		int rowsN = ceil((maxCopyNumber+1.0)/2.0);
		this->cnGenotypeBAF[cn] = ArrayXd(rowsN);
		for(int mNum=0; mNum<(cn + 2)/2; mNum++){
			double muT = 0.0;
			string piT = "";

			int pNum = cn - mNum;
			if(0 == pNum && 0 == mNum){
				muT = EMPIRI_BAF/(EMPIRI_AAF + EMPIRI_BAF);
				piT = "NULL";
			}else if(pNum == mNum){
				muT = 0.5;
				piT = this->genotype(pNum, mNum);
			}else{
				muT = (mNum * 1.0) / cn;
				piT = this->genotype(pNum, mNum) +
					'/' + this->genotype(mNum, pNum);
			}
			this->cnGenotype[cn].push_back(piT);
			this->cnGenotypeBAF[cn](mNum) = muT;
		}
	}
}
