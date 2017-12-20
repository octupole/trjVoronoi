#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([trjVoronoi],[0.2],[Massimo.Marchi@cea.fr])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/trjVoronoi.cpp])

AC_CONFIG_MACRO_DIR([m4])

AC_PREFIX_DEFAULT($HOME/trjlibs)

AC_PROG_RANLIB
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CXX_C_O 
AC_PROG_CXXCPP
AC_PROG_INSTALL
# check for maximal optimization 
AX_COMPILER_VENDOR()
AX_CXX_MAXOPT()

# Checks for libraries.



AC_ARG_WITH(xdr-dir,
        AS_HELP_STRING([--with-xdr-dir=dir],[Provide the xdrfile root directory [xdrfile-src/xdr] ]),
        [xdr_dir=$withval],[xdr_dir="$PWD/xdrfile-src/xdr"])


CPPFLAGS="$CPPFLAGS -I$xdr_dir/include -I../trjlib -I../Molecules -I../voro++ -I../Voronoi -I../Properties -I../alglib "
LDFLAGS="$LDFLAGS  -L$xdr_dir/lib -L../Voronoi -L../voro++ -L../Molecules -L../trjlib -L../Properties -L../alglib"

AC_SEARCH_LIBS([read_xtc_natoms],[xdrfile xdrfile_mpi],[],AC_MSG_ERROR([xdrfile is needed. You must provide a path to the libraryor use --with-xdr-dir=dir]))

# Checks for header files.

AC_HEADER_STDC   

AC_CHECK_HEADERS([xdrfile.h],,AC_MSG_ERROR([No header files found for XDRFILE! You must provide a path to the include files or the root installation directory with --with-xdr-dir=root]))


# Check whether the user wants to use a particular MPI installation

AC_ARG_ENABLE(mpi, [AS_HELP_STRING([--enable-mpi],[compile trjSaxs with MPI support])], enable_mpi=$enableval, enable_mpi=no)


if test x"$enable_mpi" = x"yes"; then
   AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
   AC_CHECK_PROGS(MPICXX, mpic++ mpiCC mpicxx mpCC hcp mpxlC mpxlC_r cmpic++, "no")
   AC_SUBST(MPICXX)
   if test x"$MPICXX" = x"no"; then
      AC_MSG_ERROR([mpi is needed])
   fi
   CXX=$MPICXX
   AC_MSG_CHECKING([Searching for MPI libraries])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[int argc; char **argv; MPI_Init(&argc,&argv);]])],[AC_MSG_RESULT(yes)],[AC_MSG_ERROR(Cannot compile and link MPI code with $MPICXX)])

   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   AC_SUBST(MPIRUN)
   AC_DEFINE([HAVE_MPI],1,[Enable MPI compilation])
fi

AM_CONDITIONAL(MPI, test "$enable_mpi" = "yes")

LIBS=" -lvoro -lvoro++  -lmol -ltrj -lprop -lalg $LIBS "

# Checks for programs.


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.

AC_CHECK_FUNCS([remove sqrt strdup mkstemp exit])

AC_CONFIG_FILES([Makefile Voronoi/Makefile Molecules/Makefile voro++/Makefile 
     trjlib/Makefile src/Makefile Properties/Makefile] alglib/Makefile)
AC_OUTPUT

