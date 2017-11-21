#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([],[0.1],[Massimo.Marchi@cea.fr])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/xdrfile_c_test.c])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_PREFIX_DEFAULT(../xdr)

AC_PROG_RANLIB
AC_LANG([C])
AC_PROG_CC
AC_PROG_CC_C_O
AC_PROG_CPP
AC_PROG_INSTALL
# check for maximal optimization 
AX_COMPILER_VENDOR()
AX_CC_MAXOPT()

# check for bool (must come late, since CFLAGS/includes might affect it)
AC_CHECK_TYPES([bool])

# Checks for libraries.




# Check whether the user wants to use a particular MPI installation

AC_ARG_ENABLE(mpi, [AS_HELP_STRING([--enable-mpi],[compile dftChi with MPI support])], enable_mpi=$enableval, enable_mpi=no)


if test x"$enable_mpi" = x"yes"; then
   AC_ARG_VAR(MPICC,[MPI C compiler command])
   AC_CHECK_PROGS(MPICC, mpicc mpicxx mpcc, "no")
   AC_SUBST(MPICC)
   if test x"$MPICC" = x"no"; then
      AC_MSG_ERROR([mpi is needed])
   fi
   CC=$MPICC
   AC_MSG_CHECKING([Searching for MPI libraries])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[int argc; char **argv; MPI_Init(&argc,&argv);]])],[AC_MSG_RESULT(yes)],[AC_MSG_ERROR(Cannot compile and link MPI code with $MPICXX)])

   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   AC_SUBST(MPIRUN)
   AC_DEFINE([__MPI],1,[Enable MPI compilation])
fi

AM_CONDITIONAL(MPI, test "$enable_mpi" = "yes")

AH_TEMPLATE([F77_OR_C_FUNC],
            [Set to F77_FUNC(name,NAME) if Fortran used, otherwise 'name' for C.])
AH_TEMPLATE([F77_OR_C_FUNC_],
            [Set to F77_FUNC_(name,NAME) if Fortran used, otherwise 'name' for C.])

  AC_DEFINE([F77_OR_C_FUNC(name,NAME)],[name])
  AC_DEFINE([F77_OR_C_FUNC_(name,NAME)],[name])

# Checks for programs.


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_FUNC_MEMCMP
AC_TYPE_SIGNAL
AC_FUNC_VPRINTF
AM_WITH_DMALLOC
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(strdup)

######
# check for xtc headers
case "${host_os}" in
cygwin* | mingw*)   # use our own xdr stuff for cygwin/mingw
  GMXLIB_LIBOBJS="$GMXLIB_LIBOBJS gmx_system_xdr.lo"
  ;;
*)
  AC_CHECK_HEADERS([rpc/rpc.h])
  AC_CHECK_HEADERS([rpc/xdr.h],,,
[#if HAVE_RPC_RPC_H
# include <rpc/rpc.h>
#endif
])
  AC_CHECK_LIB(nsl,xdr_float)
  AC_TRY_LINK([
#if HAVE_RPC_RPC_H
#include<rpc/rpc.h> 
#else
#error NO rpc.h header
#endif
#if HAVE_RPC_XDR_H
#include<rpc/xdr.h>
#else
#error No xdr.h header
#endif
],[ XDR *xd; float f; xdr_float(xd,&f);],[have_xdr="yes"])

  if test "$have_xdr" != "yes"; then
    AC_DEFINE(GMX_INTERNAL_XDR,,[Use our own instead of system XDR libraries])	
      GMXLIB_LIBOBJS="${GMXLIB_LIBOBJS} gmx_system_xdr.lo"
    AC_MSG_WARN([Couldn't find XDR headers and/or libraries - using our own])
  fi
  ;;
esac

############################################################
# Add some debug info: Who is building, and on what machine?
############################################################
AC_DEFINE_UNQUOTED([BUILD_TIME],["`date`"],[Date and time for build])
AC_DEFINE_UNQUOTED([BUILD_USER],["`whoami`@`hostname`"],[User doing build])
AC_DEFINE_UNQUOTED([BUILD_MACHINE],["`uname -srm`"],[Hardware and OS version for build host])

#############
# Check floating-point format and endian
#############
ACX_FLOAT_FORMAT

#############
# Check integer endian
#############
AC_C_BIGENDIAN([AC_DEFINE([GMX_INTEGER_BIG_ENDIAN],,[Integer byte order is big endian.])],[],[
AC_MSG_ERROR([Cannot determine endian in compiled output])])
# Little endian implied if GMX_INTEGER_BIG_ENDIAN not defined in source.

#############
# Check integer sizes, set to zero if e.g. long long is not supported.
#############
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(long long int)


#############
# Make sure size_t can hold pointers.
# On some Nec machines, you reportedly need to specify a flag for size_t to be 64 bits wide.
#############
AC_MSG_CHECKING([that size_t can hold pointers])
AC_TRY_COMPILE([
#include <stdlib.h>
#include <stddef.h>
],[int chk[sizeof(size_t)-sizeof(int *)+1];],AC_MSG_RESULT([yes]),[
   AC_MSG_ERROR([size_t is smaller than int*. Check if your compiler needs a flag to make it 64 bit.])])


AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Checks for library functions.

AC_CHECK_FUNCS([remove sqrt strdup mkstemp exit])

AC_CONFIG_FILES([Makefile include/Makefile src/Makefile])
AC_OUTPUT

