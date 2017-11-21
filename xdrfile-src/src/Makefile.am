# Makefile for 
# Copyright (C) 2012 Massimo Marchi <Massimo.Marchi@cea.fr>
# $Header$
# $Log$
# what flags you want to pass to the C compiler & linker
AM_CPPFLAGS = -I$(top_srcdir)/include


# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = trr2xtc

lib_LIBRARIES = libxdrfile.a
libxdrfile_a_SOURCES = xdrfile.c xdrfile_trr.c xdrfile_xtc.c	xdrfile_seek.c

LDADD = libxdrfile.a
