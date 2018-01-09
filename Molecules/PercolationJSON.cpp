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
	std::hash<string> str_hash;
	for(size_t o=0;o<this->Clusters.size();o++){
		json myClust;
		mapRes.clear();mapResAtm.clear();
		std::map<string,vector<int>> Tag;
		for(size_t p=0;p<this->Clusters[o].size();p++){
			int n=this->Clusters[o][p];
			mapRes[this->pResn[n]]++;
			Tag[this->pResn[n]].push_back(n);
			for(size_t q{0};q<this->Atoms[n].size();q++){
				mapResAtm[this->pResn[n]]++;
			}
		}
		json tmp;
		for(auto it{mapRes.begin()};it!= mapRes.end();it++){
			std::stringstream ss;
			std::copy(Tag[it->first].begin(),Tag[it->first].end(),std::ostream_iterator<int>( ss," "));
			tmp[it->first]={it->second,mapResAtm[it->first],str_hash(ss.str())};
		}
		myJson["cluster"].push_back(tmp);
	}

	fout <<"{\"cluster\": ";
	fout<<myJson["cluster"]<<",";
	myJson.clear();
}

template class PercolationJSON<float>;
template class PercolationJSON<double>;
