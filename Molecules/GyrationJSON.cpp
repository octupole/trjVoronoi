/*
 * GyrationJSON.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: marchi
 */

#include "GyrationJSON.h"
template <typename T>
set<string> GyrationJSON<T>::Time;


template <typename T>
void GyrationJSON<T>::__Writeit(ostream & fout, string label, int o){
	static bool firstTime{true};
	string TimeC=std::to_string(this->time_c);

	json & myType=this->myJson[label];
	json myClust;
	myClust["Rg"]=sqrt(this->Radg);
	myClust["I"].push_back(this->I[XX]);
	myClust["I"].push_back(this->I[YY]);
	myClust["I"].push_back(this->I[ZZ]);
	myClust["ax"].push_back(this->axis[XX]);
	myClust["ax"].push_back(this->axis[YY]);
	myClust["ax"].push_back(this->axis[ZZ]);
	myType.push_back(myClust);
	if(!Time.count(TimeC)){
		Time.insert(TimeC);
		if(firstTime){
			firstTime=false;
			fout <<"{"<<endl;
			fout <<"\""<<TimeC<<"\": ";
		}else{
			fout <<"\""<<"Gyr"<<"\": ";
			fout<< this->myJson;
			this->myJson.clear();
			fout <<"},";
			fout <<"\""<<TimeC<<"\": ";
		}

	}
}

template <typename T>
GyrationJSON<T>::~GyrationJSON() {
	// TODO Auto-generated destructor stub
}

template class GyrationJSON<float>;
template class GyrationJSON<double>;
