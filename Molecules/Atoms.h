/*
 * Atoms.h
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */

#ifndef ATOMS_H_
#define ATOMS_H_

#include "Metric.h"
#include "AtomIndex.h"
#include<iostream>
#include <cmath>
#include <vector>
#include "MyUtilClass.h"
#include "HeaderTrj.h"
#include "Topol.h"
#include "Contacts.h"

#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_seek.h"
#include "Percolation.h"

using namespace DVECT;
using namespace MATRIX;
using std::vector;
using namespace Enums;


template <typename T>
class brot{
	MMatrix<T> co;
public:

	const MMatrix<T> & operator()(const double a,const double b,const double c,
			const double alfa,const double beta,const double gamma){
	double degrad=M_PI/180.0;

	double ax=a;
    double alf=cos(degrad*alfa);
    double bet=cos(degrad*beta);
    double qt=sin(degrad*gamma);
    double gam=cos(degrad*gamma);
    double bx=b*gam;
    double by=b*qt;
    double cx=c*bet;
    double cy=c*(alf-bet*gam)/qt;
    double cz=sqrt(c*c-cx*cx-cy*cy);
    co[YY][XX]=0.0;
    co[ZZ][XX]=0.0;
    co[ZZ][YY]=0.0;
    co[XX][XX]=ax;
    co[XX][YY]=bx;
    co[XX][ZZ]=cx;
    co[YY][YY]=by;
    co[YY][ZZ]=cy;
    co[ZZ][ZZ]=cz;
    return co;
	}

};

template <typename T>
class Atoms {
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
protected:
	bool firsttime{true};
	int nr{0};
	int status{0};
	vector<Dvect> x;
	vector<Dvect> xa;
	Metric<T> Mt;
	static int    step_c;
	static float  prec_c,time_c;
	vector<double> rd;

	virtual void ReadaStep(FstreamC * );
	virtual void ReadaStep(FstreamF * );
	virtual void ReadaStep(std::ifstream &);
	virtual void ReadaStepF(std::ifstream &);
	virtual void moveOffset(FstreamC *);
	virtual void moveOffset(FstreamF *);
	virtual void moveOffset(std::ifstream &);
	virtual void WriteaStep(FstreamC * );
	virtual void WriteaStep(FstreamF * );
	virtual void __ReconstructOneCluster(vector<bool> & a){};
	virtual Dvect __FindCell(const vector<vector<int>> & x,const vector<vector<int>> & y){return Dvect();};
	Percolation<T> * Perco{nullptr};

public:
	Atoms(){};
	Atoms(const int);
	Atoms(const AtomIndex &);
	Atoms(const Atoms &);
	virtual ~Atoms(){delete Perco;};
	struct plane{
		Typedefs::real xc[4];
		int n;
		plane & operator=(Typedefs::rvec & xa){for(int o=0;o<DIM;o++) xc[o]=xa[o];return *this;};
		plane & operator=(Typedefs::real dd){xc[3]=dd;return *this;};
		plane & operator=(double dd){xc[3]=dd;return *this;};
		plane & operator=(int i){n=i;return *this;};
	} ax;

	void setDim(const int n);
	void setCoord(const Metric<T> &, const rvec *, const AtomIndex & );
	void setCoord(const Metric<T> &, const rvec *);
	void setMT(const Metric<T> &);
	Atoms &  operator=(const Atoms &);
	Atoms & operator=(const T);
	Dvect & operator[](const int i){return x[i];};
	Atoms & operator()(const int);

	void Rot(const Matrix);
	int getNR()const{return nr;};
	const Metric<T> & getMt() const {return Mt;};
	Atoms COtoOC();
	Atoms OCtoCO();
	void doCOtoOC();
	void doOCtoCO();
	Atoms Shift(const rvec);
	Typedefs::real Dist(const int,const int);
	vector<Dvect> getX() const {return x;};
	vector<Dvect> getXA() const {return xa;};
	void setrd(vector<double> & rdx) {rd=rdx;};
	void setrd(Topol_NS::Topol & y){
		rd=y.getrd();

	};
	double getrd(int n){return rd[n];}

	float getTime(){return time_c;}
	void SetupPercolate();
	void SetupPercolate(Topol_NS::Topol &x);
	int Percolate();
	int Percolate(double y){
		this->Perco->setRcut(y);
		return this->Percolate();
	}

	bool hasPerco() const {if(Perco) return true;return false;}
	Percolation<T> * gPerco(){return Perco;}

	void pdb(const vector<string> & c);

	friend Fstream & operator+=(Fstream & fin, Atoms & y){
		if(FstreamC * finC=dynamic_cast<FstreamC *> (&fin))
			y.moveOffset(finC);
		else if(FstreamF * finF=dynamic_cast<FstreamF *> (&fin))
			y.moveOffset(finF);

		return fin;
	}

	friend std::ifstream & operator+=(std::ifstream & fin,Atoms & y){
		y.moveOffset(fin);
		return fin;
	}

	friend Fstream & operator>>(Fstream & fin, Atoms & y){
		if(FstreamC * finC=dynamic_cast<FstreamC *> (&fin))
			y.ReadaStep(finC);
		else if(FstreamF * finF=dynamic_cast<FstreamF *> (&fin))
			y.ReadaStep(finF);
		return fin;
	}
	friend Fstream & operator<<(Fstream & fout, Atoms & y){
		if(FstreamC * foutC=dynamic_cast<FstreamC *> (&fout))
			y.WriteaStep(foutC);
		return fout;
	}

	friend std::ifstream & operator>>(std::ifstream & fin, Atoms & y){
		y.ReadaStep(fin);
		return fin;
	}

};

#endif /* ATOMS_H_ */

