# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -std=c++11 



# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libprop.a 
libprop_a_SOURCES = CenterMass.cpp CenterMassMicelles.cpp CenterMassBW3.cpp CenterMassWigner.cpp\
	Quaternions.cpp QuaternionUtilities.cpp WignerDMatrices.cpp Combinatorics.cpp CMPick.cpp
	
LDADD = libprop.a
include_HEADERS = CenterMass.h CenterMassMicelles.h CenterMassBW3.h CenterMassWigner.h\
	Quaternions.hpp QuaternionUtilities.hpp WignerDMatrices.hpp Combinatorics.hpp CMPick.h
