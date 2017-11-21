# trjVoronoi
## Synopsis

trjVoronoi computes the atomic Voronoi volumes and surfaces of a molecular dynamic trajectory 
written in .dcd or .xtc format. It requires a .pdb file containing the coordinates of the whole system, including waters and ions. 
The code uses voro++ by Chris H. Rycroft to compute properties of molecular systems in solution from molecular dynamics trajectories 
obtained with GROMACS and NAMD. 
## Code Example


## Motivation

This project was created to compute the volumes of proteins residues, it then evolved to calculation of the solvent accessible surface 

## Installation
trjVoronoi uses a standard autoconf script and makefiles created by automake, like most GNU programs. 
This means your normal installation actions will be limited to run a ./configure, a make and a make install. 
The code is dependent on the following packages:

* A modified xdrfile library, which is provided with the distribution. It must be installed before the main program and it 
is contained in the directory $SOURCE/xdrfile-src, where $SOURCE is the directory where the trjSaxs source is. 


For installation, after you have changed directory to $SOURCE/xdrfile-src, issue 

sh bootstrap <br />
./configure<br />
make<br />
make install<br />

Following the installation of the xdr package, you should return to the installation directory and issue:<br />
sh bootstrap <br />
./configure <br />
The command ./configure --help will tell what are the options available for compilation. Noticeably ./configure 
has amongst others the following options:

* --enable-intel          compile without user-define literals to compile on an intel compiler<br />
* --enable-mpi           compile with MPI support<br />

Finally, the installation is completed by: <br/>
make<br />
make install<br />

trjVoronoi uses several features available with the C++11 standard, thus it will not compile with GNU g++ compilers earlier 
than 4.8 and with Intel icpc compilers versions earlier than 15.0.2. It has been successfully compiled on Centos, Debian and 
Mac OS X systems supporting those versions of the C++ compilers.

## How to cite trjVoronoi: 

Please consider first citing the Voro++  code with one of the three references:

Chris H. Rycroft, Voro++: A three-dimensional Voronoi cell library in C++, Chaos 19, 041111 (2009).
Chris H. Rycroft, Gary S. Grest, James W. Landry, and Martin Z. Bazant, Analysis of Granular Flow in a Pebble-Bed Nuclear Reactor, Phys. Rev. E 74, 021306 (2006).
Chris H. Rycroft, Multiscale Modeling in Granular Flow, PhD thesis submitted to the Massachusetts Institute of Technology, September 2007.

Then consider citing our Work:
S. Abel, F. Y. Dupradeau and M. Marchi "Molecular Dynamics Simulations of a Characteristic DPC Micelle in Water" J. Chem. Theory Comp. 8, 4610-4623 (2012)

## License

  CeCILL FREE SOFTWARE LICENSE AGREEMENT

Version 2.1 dated 2013-06-21


This Agreement is a Free Software license agreement that is the result
of discussions between its authors in order to ensure compliance with
the two main principles guiding its drafting:

  * firstly, compliance with the principles governing the distribution
    of Free Software: access to source code, broad rights granted to users,
  * secondly, the election of a governing law, French law, with which it
    is conformant, both as regards the law of torts and intellectual
    property law, and the protection that it offers to both authors and
    holders of the economic rights over software.

Find the full licence here: http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt
