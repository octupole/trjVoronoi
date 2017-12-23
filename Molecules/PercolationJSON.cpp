/*
 * PercolationJSON.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: marchi
 */

#include "PercolationJSON.h"

template <typename T>
PercolationJSON<T>::~PercolationJSON() {
	// TODO Auto-generated destructor stub
}

template <typename T>
void PercolationJSON<T>::__Writeit(ostream & fout){
	map<string,int> mapRes;
	map<string,int> mapResAtm;
	for(size_t o=0;o<this->Clusters.size();o++){
		json myClust;
		mapRes.clear();mapResAtm.clear();
		for(size_t p=0;p<this->Clusters[o].size();p++){
			int n=this->Clusters[o][p];
			mapRes[this->pResn[n]]++;
			for(size_t q{0};q<this->Atoms[n].size();q++){
				mapResAtm[this->pResn[n]]++;
			}
		}
		json tmp;
		for(auto it{mapRes.begin()};it!= mapRes.end();it++){
			tmp[it->first]={it->second,mapResAtm[it->first]};
		}
		myJson["cluster"].push_back(tmp);
	}

	fout <<"{\"cluster\": ";
	fout<<myJson["cluster"]<<",";
	myJson.clear();
}

template class PercolationJSON<float>;
template class PercolationJSON<double>;
