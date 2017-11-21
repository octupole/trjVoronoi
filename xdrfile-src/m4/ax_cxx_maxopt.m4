dnl @synopsis AX_CC_MAXOPT
dnl @summary turn on optimization flags for the C compiler
dnl @category C
dnl
dnl Try to turn on "good" C optimization flags for various compilers
dnl and architectures, for some definition of "good".  (In our case,
dnl good for FFTW and hopefully for other scientific codes.  Modify 
dnl as needed.)
dnl
dnl The user can override the flags by setting the CFLAGS environment
dnl variable.  
dnl
dnl Note also that the flags assume that ANSI C aliasing rules are
dnl followed by the code (e.g. for gcc's -fstrict-aliasing), and that
dnl floating-point computations can be re-ordered as needed.
dnl
dnl Requires macros: AX_CHECK_COMPILER_FLAGS, AX_COMPILER_VENDOR,
dnl
dnl @version 2011-06-22
dnl @license GPLWithACException
dnl @author Steven G. Johnson <stevenj@alum.mit.edu> and Matteo Frigo.
AC_DEFUN([AX_CXX_MAXOPT],
[
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])
# Try to determine "good" native compiler flags if none specified via CFLAGS
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS=""
  case $ax_cv_cxx_compiler_vendor in
    dec) CXXFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
    	 ;;

    sun) CXXFLAGS="-native -fast -xO5 -dalign"
    	 ;;

    hp)  CXXFLAGS="+Oall +Optrs_ansi +DSnative"
    	 ;;

    ibm) xlc_opt="-qtune=auto"
         AX_CHECK_COMPILER_FLAGS($xlc_opt,
         	CXXFLAGS="-O3 -qansialias -w $xlc_opt",
               [CXXFLAGS="-O3 -qansialias -w"
                echo "******************************************************"
                echo "*  You seem to have the IBM  C compiler.  It is      *"
                echo "*  recommended for best performance that you use:    *"
                echo "*                                                    *"
                echo "*    CXXFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
                echo "*                      ^^^        ^^^                *"
                echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
                echo "*  CPU you have.  (Set the CXXFLAGS environment var.   *"
                echo "*  and re-run configure.)  For more info, man cc.    *"
                echo "******************************************************"])
         ;;

    intel) CXXFLAGS="-O3"
        # Intel seems to have changed the spelling of this flag recently
        icc_ansi_alias="unknown"
	for flag in -ansi-alias -ansi_alias; do
	  AX_CHECK_COMPILER_FLAGS($flag, [icc_ansi_alias=$flag; break])
	done
 	if test "x$icc_ansi_alias" != xunknown; then
            CXXFLAGS="$CXXFLAGS $icc_ansi_alias"
        fi
	AX_CHECK_COMPILER_FLAGS(-malign-double, CXXFLAGS="$CXXFLAGS -malign-double")
	# We used to check for architecture flags here, e.g. -xHost etc.,
	# but these flags are problematic.  On icc-12.0.0, "-mavx -xHost"
	# overrides -mavx with -xHost, generating SSE2 code instead of AVX
	# code.  ICC does not seem to support -mtune=host or equivalent
	# non-ABI changing flag.
	;;
    
    gnu) 
     # Default optimization flags for gcc on all systems.
     # Somehow -O3 does not imply -fomit-frame-pointer on ia32
     CXXFLAGS="-O3 -fomit-frame-pointer"

     # tune for the host by default
     AX_CHECK_COMPILER_FLAGS(-mtune=native, CXXFLAGS="$CXXFLAGS -mtune=native")

     # -malign-double for x86 systems
     AX_CHECK_COMPILER_FLAGS(-malign-double, CXXFLAGS="$CXXFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILER_FLAGS(-fstrict-aliasing,
	CXXFLAGS="$CXXFLAGS -fstrict-aliasing")

     # -fno-schedule-insns is pretty much required on all risc
     # processors.
     # 
     # gcc performs one pass of instruction scheduling, then a pass of
     # register allocation, then another pass of instruction
     # scheduling.  The first pass reorders instructions in a way that
     # is pretty much the worst possible for the purposes of register
     # allocation.  We disable the first pass.
     AX_CHECK_COMPILER_FLAGS(-fno-schedule-insns, CXXFLAGS="$CXXFLAGS -fno-schedule-insns")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILER_FLAGS(-ffast-math, CXXFLAGS="$CXXFLAGS -ffast-math")

     ;;
  esac

  if test -z "$CXXFLAGS"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best CXXFLAGS for this system  *"
        echo "* Use ./configure CXXFLAGS=... to specify your own flags *"
	echo "* (otherwise, a default of CXXFLAGS=-O3 will be used)    *"
	echo "********************************************************"
	echo ""
        CXXFLAGS="-O3"
  fi

  AX_CHECK_COMPILER_FLAGS($CXXFLAGS, [], [
	echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CXXFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure CXXFLAGS=... to specify your own flags *"
        echo "********************************************************"
        echo ""
        CXXFLAGS=""
  ])

fi
])
