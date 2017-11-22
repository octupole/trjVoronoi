# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -std=c++11 



# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libtrj.a
libtrj_a_SOURCES = FComms.cpp  Fstream.cpp  FstreamF.cpp  Metric.cpp  NewMPI.cpp  Timer.cpp\
	Finalize.cpp  FstreamC.cpp  HeaderTrj.cpp  MyUtilClass.cpp  Split.cpp  histograms.cpp

LDADD = libtrj.a
include_HEADERS = Array.h Fstream.h  HeaderTrj.h  NewMPI.h \
	Communicator.hpp FstreamC.h  Metric.h  PickSelection.h  histograms.hpp\
	FComms.h  FstreamF.h  MyUtilClass.h  Split.h   myEnums.hpp\
	Finalize.h  Ftypedefs.h  MyVec.h   Timer.h

