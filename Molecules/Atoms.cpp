/*
 * Atoms.cpp
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */
#include "Atoms.h"
#include "jacobi.h"
template <typename T>
class Jacob{
	double ** iner;
	double ** ei;
	double * Im;
	static Dvect Im_s;
	void  eigsort(DDvect<T>  & d, MMatrix<T>  * v=NULL){
		int k;
		int n=DIM;
		for (int i=0;i<n-1;i++) {
			double p=d[k=i];
			for (int j=i;j<n;j++)
				if (d[j] >= p) p=d[k=j];
			if (k != i) {
				d[k]=d[i];
				d[i]=p;
				if (v != NULL){
					for (int j=0;j<n;j++) {
						p=(*v)[j][i];
						(*v)[j][i]=(*v)[j][k];
						(*v)[j][k]=p;
					}
				}
			}
		}
	}
public:
	Jacob():  iner(NULL), ei(NULL),Im(NULL){};
	~Jacob(){
		for(int o=0;o<DIM;o++) delete [] iner[o];
		delete [] iner;
		for(int o=0;o<DIM;o++) delete [] ei[o];
		delete [] ei;
		delete [] Im;
	}
	void operator()(MMatrix<T> & iner0, DDvect<T> & Im0, MMatrix<T> & ei0, int & nrot){
		if(!iner) {
			iner=new double * [DIM];
			for(int o=0;o<DIM;o++) iner[o]=new double [DIM];
		}
		if(!ei) {
			ei=new double * [DIM];
			for(int o=0;o<DIM;o++) ei[o]=new double [DIM];
		}
		if(!Im) Im=new double [DIM];
		for(int o=0;o<DIM;o++){
			for(int p=0;p<DIM;p++){
				iner[o][p]=iner0[o][p];
			}
		}
		jacobi(iner,DIM,Im,ei,&nrot);

		for(int o=0;o<DIM;o++){
			Im0[o]=Im[o];
			for(int p=0;p<DIM;p++){
				ei0[o][p]=ei[o][p];
			}
		}

		eigsort(Im0,&ei0);
	}

};
template <typename T>
Dvect Jacob<T>::Im_s=0.0;

template<typename T>
int    Atoms<T>::step_c=0;

template<typename T>
float   Atoms<T>::prec_c=0.0;

template<typename T>
float Atoms<T>::time_c=0.0;

template <typename T>
void Atoms<T>::setTopol(Topol_NS::Topol & y){
	const double massH{1.00784};
	const double massC{12.00960};
	const double eps{1.0e-4};
	rd=y.getrd();
	mass=y.getMass();
	massNCH=mass;
	std::transform(massNCH.begin(),massNCH.end(),massNCH.begin(),[=](double x){
		return (abs(x-massH) < eps || abs(x-massC) <eps)?0.0:x;
	});
	for(size_t o{0};o<mass.size();o++){

	}
};

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
void Atoms<T>::CalcGyro(vector<double> & massa,vector<Gyration<T>> & Rg){
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();
	for(size_t o=0;o<mCluster.size();o++){
		Rg[o]=0.0;
		Dvect cm,cmG;
		double unit_nmm2=1.0/(unit_nm*unit_nm);

		int ntot=0;
		double tmass=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int i=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) cm[o1]+=x[i][o1]*massa[i];
				for(int o1=0;o1<DIM;o1++) cmG[o1]+=x[i][o1];
				tmass+=massa[i];
				ntot++;
			}
		}
		cm/=tmass;
		cmG/=static_cast<int> (ntot);
		Dvect x0;
		Dvect Gm,Im,axis;
		Matrix Giner;
		Matrix ei,ei2;
		double MyRg=0;
		int oa=0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int nn=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) x0[o1]=cmG[o1]-x[nn][o1];
				Giner+=x0%x0*unit_nmm2;
			}
		}
		Giner/=static_cast<double> (ntot);
		int nrot;
		Jacob<T>()(Giner,Gm,ei,nrot);
		for(int o1=0;o1<DIM;o1++) MyRg+=Gm[o1];
		Matrix Inertia;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int nn=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) x0[o1]=cm[o1]-x[nn][o1];

				Inertia[XX][XX]+=(sqr(x0[YY])+sqr(x0[ZZ]))*massa[nn];
				Inertia[YY][YY]+=(sqr(x0[ZZ])+sqr(x0[XX]))*massa[nn];
				Inertia[ZZ][ZZ]+=(sqr(x0[XX])+sqr(x0[YY]))*massa[nn];
				Inertia[XX][YY]-=(x0[XX]*x0[YY])*massa[nn];
				Inertia[XX][ZZ]-=(x0[XX]*x0[ZZ])*massa[nn];
				Inertia[YY][ZZ]-=(x0[YY]*x0[ZZ])*massa[nn];
			}
		}
		Inertia[YY][XX]=Inertia[XX][YY];
		Inertia[ZZ][XX]=Inertia[XX][ZZ];
		Inertia[ZZ][YY]=Inertia[YY][ZZ];
		Jacob<T>()(Inertia,Im,ei2,nrot);
		Im*=unit_nmm2;
		double fact=5.0/tmass/2.0;
		axis[XX]=fact*(Im[YY]+Im[ZZ]-Im[XX]);
		axis[YY]=fact*(Im[XX]+Im[ZZ]-Im[YY]);
		axis[ZZ]=fact*(Im[YY]+Im[XX]-Im[ZZ]);
		MyRg=(axis[XX]+axis[YY]+axis[ZZ])/5.0;
		Rg[o](MyRg,Im,Gm,axis);
	}

}
template <typename T>
void Atoms<T>::Gyro(){
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();

	vector<Gyration<T>> Rg=vector<Gyration<T>>(mCluster.size());
	Rg_i.clear();
	Gyration<T>::setTime(time_c);
	CalcGyro(mass,Rg);
	copy(Rg.begin(),Rg.end(),back_inserter(Rg_i));
	CalcGyro(massNCH,Rg);
	copy(Rg.begin(),Rg.end(),back_inserter(Rg_i));
	Rg_count++;
}

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
void Atoms<T>::SetupPercolate(bool JSON){
}
template <typename T>
void Atoms<T>::SetupPercolate(Topol_NS::Topol & myTop, bool JSON){
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
template <typename T>
vector<DDvect<T>> Atoms<T>::getGC(){
	try{
		if(!this->hasPerco()) throw string("Cannot compute clusters center of mass without percolation.");
	}catch(const string & s){
		cout << s<<endl;exit(1);
	}
	const vector<vector<int> > & mCluster=Perco->getCluster();
	vector<vector<int> > & mAtoms=Perco->getAtoms();
	vector<Dvect> R_CM=vector<Dvect>(mCluster.size());
	for(size_t o=0;o<mCluster.size();o++){
		Dvect cm{0};
		T tmass=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int i=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) cm[o1]+=xa[i][o1];
				tmass+=1.0;
			}
		}
		cm/=tmass;
		R_CM[o]=cm;
	}
	return R_CM;
}

template class Atoms<float>;
template class Atoms<double>;

