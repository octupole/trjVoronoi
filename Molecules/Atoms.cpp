/*
 * Atoms.cpp
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */
#include "Atoms.h"

template<typename T>
int    Atoms<T>::step_c=0;

template<typename T>
float   Atoms<T>::prec_c=0.0;

template<typename T>
float Atoms<T>::time_c=0.0;


template <typename T>
Atoms<T>::Atoms(const int nind) {
	nr=nind;
	if(nr) {
		x=vector<Dvect>(nr);
		xa=vector<Dvect>(nr);
	}
}
template <typename T>
Atoms<T>::Atoms(const AtomIndex & id){
	nr=id.getN();
	if(nr) {
		x=vector<Dvect>(nr);
		xa=vector<Dvect>(nr);
	}
}
template <typename T>
Atoms<T>::Atoms(const Atoms<T> & y){
	nr=y.nr;
	x=vector<Dvect>(nr);
	xa=vector<Dvect>(nr);
	for(int i=0;i<nr;i++) {
		x[i][0]=y.x[i][0];
		x[i][1]=y.x[i][1];
		x[i][2]=y.x[i][2];
		xa[i][0]=y.xa[i][0];
		xa[i][1]=y.xa[i][1];
		xa[i][2]=y.xa[i][2];
	}
	Mt(y.Mt);
}

template <typename T>
Atoms<T> & Atoms<T>::operator()(const int natoms){
	nr=natoms;
	if(!x.empty()) x.clear();
	if(!xa.empty()) xa.clear();
	x=vector<Dvect>(nr);
	xa=vector<Dvect>(nr);

	return *this;
};

template <typename T>
Atoms<T> & Atoms<T>::operator=(const Atoms<T> & y){
	try{
		if(nr != y.nr) throw "Cannot copy Atoms of different size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(int i=0;i<nr;i++) {
		x[i][0]=y.x[i][0];
		x[i][1]=y.x[i][1];
		x[i][2]=y.x[i][2];
	}
	Mt(y.Mt);
	return *this;
}
template <typename T>
Atoms<T> & Atoms<T>::operator=(const T a){
	for(int i=0;i<nr;i++) {
		x[i][0]=a;
		x[i][1]=a;
		x[i][2]=a;
	}
	return *this;
}
template <typename T>
void Atoms<T>::Rot(const Matrix coR){
	rvec t;
	for(int i=0;i<nr;i++){
		for(int o=0;o<DIM;o++) t[o]=coR[o][0]*x[i][0]+coR[o][1]*x[i][1]+coR[o][2]*x[i][2];
		for(int o=0;o<DIM;o++) x[i][o]=t[o];
	}
}
template <typename T>
void Atoms<T>::setCoord(const Metric<T> & Mt_in, const rvec * x0){
	Mt(Mt_in);
	for(int i=0;i<nr;i++)
		for(int j=0;j<DIM;j++)
			x[i][j]=x0[i][j];
}
template <typename T>
void Atoms<T>::setMT(const Metric<T> & Mt_in){
		Mt(Mt_in);
	}

template <typename T>
void Atoms<T>::doCOtoOC(){
	try{
	if(x.empty()) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	if(xa.empty()) xa=vector<Dvect>(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) xa[i][o]=Mt.getOC()[o][0]*x[i][0]+Mt.getOC()[o][1]*x[i][1]+
			Mt.getOC()[o][2]*x[i][2];
}
template <typename T>
void Atoms<T>::doOCtoCO(){
	try{
	if(xa.empty() || x.empty()) throw "Atoms class atom coordinates undefined";}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) x[i][o]=Mt.getCO()[o][0]*xa[i][0]+Mt.getCO()[o][1]*xa[i][1]+
			Mt.getCO()[o][2]*xa[i][2];
}
template <typename T>
Atoms<T> Atoms<T>::COtoOC(){
	try{
	if(x.empty()) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	Atoms<T> other(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) other.x[i][o]=Mt.getOC()[o][0]*x[i][0]+Mt.getOC()[o][1]*x[i][1]+Mt.getOC()[o][2]*x[i][2];
	return other;
}
template <typename T>
Atoms<T> Atoms<T>::OCtoCO(){
	try{
	if(x.empty()) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	Atoms<T> other(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) other.x[i][o]=Mt.getCO()[o][0]*x[i][0]+Mt.getCO()[o][1]*x[i][1]+Mt.getCO()[o][2]*x[i][2];
	return other;
}
template <typename T>
void Atoms<T>::WriteaStep(FstreamF * fout){
	try{
		throw string("Don't know how to write dcd files. Maybe I should take a class... Abort");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
}
template <typename T>
void Atoms<T>::WriteaStep(FstreamC * fout){
	XDRFILE * xd=fout->getfin();


	float prec=prec_c;
	matrix box;
	rvec * x0=new rvec[nr];
	for(int o=0;o<nr;o++)
		for(int q=0;q<DIM;q++)
			x0[o][q]=static_cast<float>(x[o][q]);


	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++)
			box[i][j]=Mt.getCO()[i][j];
	try{
		if(write_xtc(xd,nr,step_c,time_c,box,x0,prec)) throw string("Cannot write next frame. Abort");
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}
	delete [] x0;

}

template <typename T>
void Atoms<T>::ReadaStep(std::ifstream & fin){
	std::cout << " Why am I here?? "<<  std::endl;
	exit(1);
}

template <typename T>
void Atoms<T>::ReadaStep(FstreamC * fin){
	XDRFILE * xd=fin->getfin();
	rvec   *x0;
	matrix box;
	int    step;
	float   prec,time;
	int my_nframe=fin->gFrame();

	FILE * fp=xdrfile_get_fp(xd);
	XDR * myxdr=xdrfile_get_xdr(xd);

	try{
		if(firsttime){
			firsttime=false;
			if(my_nframe != fin->goffStep()){
				if(!(my_nframe<= fin->gFrameNumber()))
					throw string(" Beginning frame is out of range ");
				}
			fin->Rewind();
		}
		x0=new rvec[nr];
		// Apparenly xdr_xtc_seek_frame does not work for the beginning frame!!
		if(my_nframe!=fin->goffStep()) xdr_xtc_seek_frame(my_nframe,fp,myxdr,nr);
		if(read_xtc(xd,nr,&step,&time,box,x0,&prec))
			throw string("Cannot read next frame!");
	}
	catch(const string & s){
		cout << s << endl;
		exit(-1);
	}

	prec_c=prec;
	time_c=time;
	step_c=step;
	MMatrix<T> box1;
	for(auto o=0;o<DIM;o++)
		for(auto p=0;p<DIM;p++)
			box1[o][p]=static_cast<T> (box[o][p]);
	Metric<T> Met(box1);
	setCoord(Met,x0);
	doCOtoOC();
	fin->nextFrame();
	delete []x0;
}

template <typename T>
void Atoms<T>::ReadaStep(FstreamF * fin){
	ReadaStepF(fin->getfin());
}

template <typename T>
void Atoms<T>::pdb(const vector<string> & data){

	vector<Dvect> xc;
	double aa=-1.0,bb=-1.0,cc=-1.0,alpha=90.0,beta=90.0,gamma=90.0;
	for(size_t i=0;i<data.size();i++){
		if(data[i].find("CRYST1") == 0){

			istringstream ss(data[i].substr(6,49));
			ss>>aa>>bb>>cc>>alpha>>beta>>gamma;
			aa*=unit_nm;
			bb*=unit_nm;
			cc*=unit_nm;
		}
	}
	try{
	if(aa < 0 || bb < 0 || cc <0) throw string("\nWarning: PDB file does not have CRYST1 keyword. ")+
			string("If you are doing calculations\n     on this structure, box is not set.\n");
	} catch(const string & s){cout << s<<endl;};
	Metric<T> Met(brot<T>()(aa,bb,cc,alpha,beta,gamma));
	for(size_t i=0;i<data.size();i++){
		if(data[i].find("ATOM") == 0 || data[i].find("HETATM") == 0 ) {
			Dvect tmp;
			std::stringstream ss(data[i].substr(30,24));
			ss>> tmp;
			xc.push_back(tmp);
		}
	}
	try{
		if(this->nr != int(xc.size())) throw string("Something wrong: No. atoms do not match");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}

	rvec * x0=new rvec[this->nr];
	for(int n=0;n<this->nr;n++){
		for(int o=0;o<DIM;o++) x0[n][o]=xc[n][o]*unit_nm;
	}
	this->setCoord(Met,x0);
	this->doCOtoOC();
	delete [] x0;
}


template <typename T>
void Atoms<T>::ReadaStepF(std::ifstream & fin){
	double table[6];
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.read(reinterpret_cast<char *> (&table),sizeof(table));
	fin.seekg(FORTRANBYTES,ios::cur);
	table[0]*=unit_nm;
	table[2]*=unit_nm;
	table[5]*=unit_nm;
	double alpha = 90.0 - asin(table[1]) * 90.0 / M_PI_2;
	double beta = 90.0 - asin(table[3]) * 90.0 / M_PI_2;
	double gamma = 90.0 - asin(table[4]) * 90.0 / M_PI_2;

	Metric<T> Met(brot<T>()(table[0],table[2],table[5],alpha,beta,gamma));
	vector<float> X[DIM];

	for(int i=0;i<DIM;i++){
		X[i].resize(nr,0.0);
		fin.seekg(FORTRANBYTES,ios::cur);
		fin.read(reinterpret_cast<char *> (&X[i][0]),sizeof(X[0][0])*nr);
		fin.seekg(FORTRANBYTES,ios::cur);
	}
	rvec * x=new rvec[nr];
	for(int n=0;n<nr;n++){
		for(int o=0;o<DIM;o++) x[n][o]=X[o][n]*unit_nm;
	}

	setCoord(Met,x);
	doCOtoOC();
}
template <typename T>
void Atoms<T>::moveOffset(FstreamC * fin){
	fin->nextFrame();
}
template <typename T>
void Atoms<T>::moveOffset(FstreamF * fin){
	moveOffset(fin->getfin());
}
template <typename T>
void Atoms<T>::moveOffset(std::ifstream & fin){
	const int DOUBLE=8;
	const int FLOAT=4;
	ios::streamoff OFFSET=FORTRANBYTES*2+6*DOUBLE+(FORTRANBYTES*2+nr*FLOAT)*DIM;
	fin.seekg(OFFSET,ios::cur);
}

template <typename T>
void Atoms<T>::SetupPercolate(){
}
template <typename T>
void Atoms<T>::SetupPercolate(Topol_NS::Topol & myTop){
	auto & Reference=myTop.gReferenceResidues();

	auto & atmss=myTop.getAtomName();
	auto & Index=myTop.gCIndex();

	vector<string> resn(nr);
	for(int o{0};o<nr;o++)
		resn[o]=myTop.AtomResidue(o);
	vector<vector<int>> MySel;
	for(size_t o{0};o<Reference.size();o++)
		MySel.push_back(Index[Reference[o]]);

	this->Perco=new Percolation<T>(MySel,rd,resn,atmss);
}


template <typename T>
int Atoms<T>::Percolate() {
	try{
		if(!this->Perco) throw string("Should initialize percolation. Abort.");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	vector<Dvect> v{x};
	Matrix co(Mt.getCO());
	Matrix oc(Mt.getOC());
	this->Perco->doContacts(v,co,oc);
	int result=this->Perco->gCluster();
	this->Perco->Accumulate();
	return result;
}


template class Atoms<float>;
template class Atoms<double>;

