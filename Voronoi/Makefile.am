# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -std=c++11 



# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libvoro.a
libvoro_a_SOURCES = Voronoi.cpp  VoronoiMicelles.cpp VoronoiSetter.cpp\
	VoronoiBinary.cpp IteratorVoronoi.cpp

LDADD = libvoro.a
include_HEADERS = Voronoi.h  VoronoiMicelles.h VoronoiSetter.h\
	VoronoiBinary.h IteratorVoronoi.h
