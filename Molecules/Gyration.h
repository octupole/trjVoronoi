/*
 * Gyration.h
 *
 *  Created on: Aug 11, 2013
 *      Author: marchi
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include "MyUtilClass.h"

#ifndef GYRATION_H_
#define GYRATION_H_
using namespace std;
using namespace DVECT;
using namespace MATRIX;

template <typename T>
class Gyration{
protected:
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;

	double Radg{0};
	Dvect I,G,axis;
	static double time_c;
	virtual void __Writeit(ofstream &,string, int);
public:
	Gyration(){};
	Gyration(double a,Dvect & b, Dvect & c, Dvect & d): Radg(a),I(b),G(c), axis(d) {};
	virtual ~Gyration();
	void operator()(double a,Dvect & b, Dvect & c, Dvect & d){
		Radg=a;I=b;G=c;axis=d;
	}
	void operator()(double a,Dvect & b, Dvect & c, Dvect & d, double tt){
		Radg=a;I=b;G=c;axis=d;time_c=tt;
	}
	static void setTime(double tt){time_c=tt;};
	Gyration operator/(int n){
		Gyration tmp;
		tmp.Radg=Radg/static_cast<double>(n);
		tmp.I=I/=static_cast<double>(n);
		tmp.G=G/static_cast<double>(n);
		tmp.axis=axis/static_cast<double>(n);
		return tmp;
	}
	Dvect Axis(){Dvect Ax=axis;for(int o=0;o<DIM;o++) Ax[o]=sqrt(Ax[o]);return Ax;}
	double Volume(){return (4.0/3.0)*M_PI*sqrt(axis[XX]*axis[YY]*axis[ZZ]);}
	void operator=(double y){
		Radg=y;
		I=y;G=y;axis=0.0;
	}
	Gyration & operator+=(Gyration & y){
		Radg+=y.Radg;
		I+=y.I;G+=y.G;axis+=y.axis;
		return *this;
	}
	template <typename G>
	friend ofstream & operator<<(ofstream & , Gyration<G> &);
	template <typename G>
	friend ofstream & operator<<(ofstream & , vector<Gyration<G>> &);
};
#endif /* GYRATION_H_ */
