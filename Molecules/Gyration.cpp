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
Gyration<T>::~Gyration(){};

template <typename T>
ofstream & operator<<(ofstream & fout , Gyration<T> & Rgs){
	fout << "  time = "<< fixed << setw(8) << setprecision(4) << Rgs.time_c;
	fout << "  Rg = "<< fixed << setw(8) << setprecision(4) << sqrt(Rgs.Radg);
	fout << "  a = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs.axis[XX]) ;
	fout << "  b = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs.axis[YY]) ;
	fout << "  c = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs.axis[ZZ]) ;
	fout << "  I_x = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs.I[XX]) ;
	fout << "  I_y = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs.I[YY]) ;
	fout << "  I_z = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs.I[ZZ]) ;
	fout << endl;
	return fout;
}

template <typename T>
ofstream & operator<<(ofstream & fout , vector<Gyration<T>> & Rgs){
	size_t beg0=0,beg1=Rgs.size()/2;
	size_t end0=Rgs.size()/2,end1=Rgs.size();
	for(size_t o=beg0;o<end0;o++){
		fout << "  ClusMIC = "<< fixed << setw(3) << o;
		fout << "  time = "<< fixed << setw(8) << setprecision(4) << Rgs[o].time_c;
		fout << "  Rg = "<< fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].Radg);
		fout << "  a = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].axis[XX]) ;
		fout << "  b = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].axis[YY]) ;
		fout << "  c = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].axis[ZZ]) ;
		fout << "  I_x = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].I[XX]) ;
		fout << "  I_y = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].I[YY]) ;
		fout << "  I_z = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].I[ZZ]) ;
		fout << endl;
	}
	int oo=0;
	for(size_t o=beg1;o<end1;o++){
		fout << "  ClusPOL = "<< fixed << setw(3) << oo++;
		fout << "  time = "<< fixed << setw(8) << setprecision(4) << Rgs[o].time_c;
		fout << "  Rg = "<< fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].Radg);
		fout << "  a = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].axis[XX]) ;
		fout << "  b = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].axis[YY]) ;
		fout << "  c = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].axis[ZZ]) ;
		fout << "  I_x = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].I[XX]) ;
		fout << "  I_y = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].I[YY]) ;
		fout << "  I_z = " << fixed << setw(8) << setprecision(4) << sqrt(Rgs[o].I[ZZ]) ;
		fout << endl;
	}
	return fout;
}
template class Gyration<float>;
template class Gyration<double>;
