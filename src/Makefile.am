# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -std=c++11 



# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile

bin_PROGRAMS = trjVoronoi trjProp

trjVoronoi_SOURCES = ExecuteVoronoi.cpp TrjRead.cpp  trjInput.cpp  trjVoronoi.cpp ClearUsage.cpp
trjProp_SOURCES = ExecuteProp.cpp TrjRead.cpp  trjInput.cpp  trjProp.cpp ClearUsage.cpp

trjVoronoi_DEPENDENCIES = ../trjlib/libtrj.a ../Molecules/libmol.a ../Voronoi/libvoro.a ../voro++/libvoro++.a\
	../Properties/libprop.a ../alglib/libalg.a

