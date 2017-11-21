# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -std=c++11  



# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libvoro++.a
libvoro___a_SOURCES =cell.cc common.cc container.cc container_prd.cc c_loops.cc pre_container.cc\
	unitcell.cc v_compute.cc v_base.cc wall.cc


LDADD = libvoro++.a
include_HEADERS = c_loops.hh  config.hh  pre_container.hh v_base.hh  wall.hh cell.hh container.hh\
	rad_option.hh  v_compute.hh  worklist.hh common.hh  container_prd.hh unitcell.hh  voro++.hh


