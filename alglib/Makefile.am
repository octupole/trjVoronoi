# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -g0 -std=c++11 


# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libalg.a
libalg_a_SOURCES = alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp fasttransforms.cpp\
	integration.cpp interpolation.cpp linalg.cpp optimization.cpp solvers.cpp specialfunctions.cpp\
	statistics.cpp

LDADD = libalg.a

include_HEADERS = alglibinternal.h alglibmisc.h ap.h dataanalysis.h diffequations.h fasttransforms.h\
	integrate_function_adapt_simpson.h integrate_function_adapt_simpson_abstract.h integration.h\
	interpolation.h linalg.h optimization.h solvers.h specialfunctions.h statistics.h stdafx.h

