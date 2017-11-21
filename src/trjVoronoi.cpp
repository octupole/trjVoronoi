//============================================================================
// Name        : trjVoronoi.cpp
// Author      : Massimo Marchi
// Version     :
// Copyright   : CeCILL copyright
// Description : Hello World in C++, Ansi-style
//============================================================================
/*! \mainpage trjSaxs
 *
 * \section intro_sec Overview
 *	trjSaxs computes the small angle X-ray scattering (SAXS) intensity of
 *	solvated globular and membrane proteins from 3D fast Fourier
 *	transforms (FFT). It combines a suitable particle meshing
 *	interpolation, similar to the one use in smooth particle mesh Ewald
 *	for electrostatics, with a uniform solvent density FFT padding scheme
 *	to obtain a convenient SAXS spectral resolution. The CPU time scaling
 *	of the method, as a function of system size is highly favorable and
 *	its application to solution of large systems such as solvated membrane
 *	proteins is undemanding. Differently from other approaches, all
 *	contributions from the simulation cell are included. This means that
 *	the subtraction of the buffer from the solution scattering intensity
 *	is straightforward and devoid of artefact due to ad hoc definitions of
 *	proximal and distal solvent intensity contributions.
 *
 *	For the time being, the provided documentation is not comprehensive, but
 *	it will be upgraded in the future.
 *
 * \section install_sec Installation
 *
 *	trjSaxs uses a standard autoconf script and makefiles created by
 *	automake, like most GNU programs. This means your normal installation
 *	actions will be limited to run a ./configure, a makeand a make
 *	install. The code is dependent on the following packages:
 *
 *	fftw3 for discrete Fourier transforms A modified version of the
 *	xdrfile library, which is provided with the distribution They must be
 *	installed on an accessible directory of your system before the package
 *	installation can begin. As for the xdrfile package, it is contained in
 *	the directory $SOURCE/xdrfile-src, where $SOURCE is the directory
 *	where the trjSaxs source is. After you have changed directory to
 *	$SOURCE/xdrfile-src, the sequence of commands:
 *
 *	./configure
 *	make
 *	make install
 *
 *	will install the xdrfile libraries and includes needed by trjSaxs in
 *	$HOME/xdr. This can be changed to a different directory with --prefix,
 *	but the installation of trjSaxs is simpler if you do not change the
 *	xdrfile root directory.  Following the installation of the xdr
 *	package, you should return to the installation directory and issue a
 *	./configure. The command ./configure --help will tell what are the
 *	options available for compilation. Noticeably ./configure has amongst
 *	others the following options:
 *
 *	--enable-intel		compile trjSaxs without user-define literals to  *				compile on an intel compiler
 *	--enable-openmp		compile trjSaxs with openmp support
 *	--enable-mpi		compile trjSaxs with MPI support
 *	--with-xdr-dir=dir	Provide the xdrfile root directory [$HOME/xdr]
 *	--with-fftw3-dir=dir	Provide the fftw3 root directory [/usr/local]
 *
 *	trjSaxs uses several features available with the C++11 standard, thus
 *	it will not compile with GNU g++ compilers earlier than 4.8 and with
 *	Intel icpc compilers versions earlier than 15.0.2. It has been
 *	installed on Centos and Debian Linux machines and on Mac OS X.
 */
#include <iostream>
#include <vector>
#include <string>

#include <iterator>
#include "TopolPDB.h"
#include "TopolMic.h"
#include "Timer.h"
#include "TrjRead.h"
#include "myEnums.hpp"
#include "Timer.h"
#include "ExecuteVoronoi.h"
#include "NewMPI.h"
#include "FAtoms.h"
#include "PickSelection.h"
#include "Finalize.h"
using namespace Topol_NS;

using namespace std;

Atoms<double> * atm{nullptr};
Voro::ExecuteVoronoi<double> * MyRun{nullptr};


TopolPDB topPDB;

#include <iostream>
using namespace std;
int main(int argc, char ** argv){
	/*
	**
	 * set communicator
	 */
	Voro::ExecuteVoronoi<double>::InitComm();// Start parallel tasks
	Finale::Finalize::CurrMPI=Voro::ExecuteVoronoi<double>::CurrCom();
	trj::TrjRead::SetComm(Voro::ExecuteVoronoi<double>::CurrCom());
	/**
	 * read input
	 */
	trj::TrjRead MyIn(argc,argv);
	MyIn.Input();
	Topol_NS::TopolMic MyTop;

	vector<string> data;
	if(MyIn.gFpdb()){
		// read pdb file to construct topology
		for(string str;getline(*MyIn.gFpdb(),str);){
			data.push_back(str);
		}
	}
	topPDB.sDetPolsegment(MyIn.gsDetResidue(),MyIn.gsPolarAtoms());
	topPDB(data);
	MyTop(topPDB,MyIn.bbIsrd());
	typedef map<string,vector<vector<int> > > mappa;
	mappa & Def=MyTop.getDef();

	MyIn.gReference()=PickSelection(MyIn.gReference()).Select<Enums::Reference>(Def,MyTop); // Pick the Reference residue

	int natoms=MyTop.Size();
	atm=new FAtoms<double>(natoms);

	MyTop.InitSelection(MyIn.gReference(),Enums::Reference);


	MyRun=new Voro::ExecuteVoronoi<double>(MyIn,MyTop);

	(*MyRun)(atm);

	return 0;
}
