# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -std=c++11 



# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libmol.a
libmol_a_SOURCES = AtomIndex.cpp  IteratorAtoms.cpp ResidueTypes.cpp\
	Atoms.cpp  Groups.cpp  LCells.cpp  Topol.cpp  TopolPDB.cpp Percolation.cpp \
	Contacts.cpp Gyration.cpp
	
LDADD = libmol.a
include_HEADERS = AtomIndex.h IteratorAtoms.h ResidueTypes.h \
	Atoms.h  Groups.h LCells.h Topol.h  TopolPDB.h Percolation.h \
	Contacts.h Gyration.h

