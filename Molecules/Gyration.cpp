/*
 * Gyration.cpp
 *
 *  Created on: Aug 11, 2013
 *      Author: marchi
 */

#include "Gyration.h"
template <typename T>
double Gyration<T>::time_c=-1.0;
template <typename T>
json Gyration<T>::myJson;

template <typename T>
Gyration<T>::~Gyration(){};

template <typename T>
void Gyration<T>::__Writeit(ostream & fout, string label, int o){
	fout << label<< fixed << setw(3) << o;
	fout << "  time = "<< fixed << setw(8) << setprecision(4) << time_c;
	fout << "  Rg = "<< fixed << setw(8) << setprecision(4) << sqrt(Radg);
	fout << "  a = " << fixed << setw(8) << setprecision(4) << sqrt(axis[XX]) ;
	fout << "  b = " << fixed << setw(8) << setprecision(4) << sqrt(axis[YY]) ;
	fout << "  c = " << fixed << setw(8) << setprecision(4) << sqrt(axis[ZZ]) ;
	fout << "  I_x = " << fixed << setw(8) << setprecision(4) << sqrt(I[XX]) ;
	fout << "  I_y = " << fixed << setw(8) << setprecision(4) << sqrt(I[YY]) ;
	fout << "  I_z = " << fixed << setw(8) << setprecision(4) << sqrt(I[ZZ]) ;
	fout << endl;

}

template <typename T>
ostream & operator<<(ostream & fout , Gyration<T> * Rgs){
	Rgs->__Writeit(fout,"Tot",0);
	return fout;
}
template <typename T>
ostream & operator<<(ostream & fout , vector<Gyration<T> *> & Rgs){
	size_t beg0=0,beg1=Rgs.size()/2;
	size_t end0=Rgs.size()/2,end1=Rgs.size();
	for(size_t o=beg0;o<end0;o++){
		Rgs[o]->__Writeit(fout,"Tot",o);
	}
	int oo=0;
	for(size_t o=beg1;o<end1;o++,oo++){
		Rgs[o]->__Writeit(fout,"Pol",oo);
	}
	return fout;
}
template class Gyration<float>;
template class Gyration<double>;
template ostream & operator<<(ostream &,vector<Gyration<float>*>&);
template ostream & operator<<(ostream &,vector<Gyration<double>*>&);
template ostream & operator<<(ostream &,Gyration<float>*);
template ostream & operator<<(ostream &,Gyration<double>*);
