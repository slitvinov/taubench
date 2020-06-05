dnl
dnl M4 Macros for batch job submission of MPI Jobs.
dnl This ain't easy folks ..
dnl


AC_DEFUN([AC_PROG_GNUPLOT],[
AC_PATH_PROG(GNUPLOT, gnuplot, gnuplot)
AC_SUBST(GNUPLOT)
])

AC_DEFUN([AC_PROG_PERL],[
AC_PATH_PROG(PERL, perl, nocommand)
if test x"$PERL" != "xnocommand"; then
    MYPERL="$PERL"
    AC_SUBST(MYPERL)
else
    AC_MSG_ERROR([perl not found in $PATH])
fi
])

AC_DEFUN([AC_PROG_SSH],[
AC_REQUIRE([AC_EXEEXT])dnl
AC_PATH_PROG(SSH, ssh$EXEEXT, nocommand)
if test "$SSH" = nocommand; then
       AC_MSG_ERROR([ssh not found in $PATH])
fi;dnl
])


AC_DEFUN([AC_PROG_RSH],[
AC_REQUIRE([AC_EXEEXT])dnl
AC_PATH_PROG(RSH, rsh$EXEEXT, nocommand)
if test "$RSH" = nocommand; then
       AC_MSG_ERROR([rsh not found in $PATH])
fi;dnl
])
dnl

AC_DEFUN([AC_MPI], [

AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
         AC_ARG_VAR(MPICC,[MPI C compiler command])
        AC_PATH_PROGS(MPICC,
		mpicc hcc mpcc mpcc_r mpxlc cmpicc sxmpicc, nocommand)
	if test x"$MPICC" = "xnocommand"; then
     	     AC_MSG_ERROR([Can't find MPICC])	
	fi
	CC="$MPICC"	
	AC_PROG_CC()
	mydir=`pwd`
	mympicc="$MPICC"
	MPI_BIN_PATH=`dirname "$MPICC"`
	cd $MPI_BIN_PATH
        while test -h "$MPICC"; do
             MPICC=`readlink "$MPICC"`
             MPI_BIN_PATH=`dirname "$MPICC"`
	     cd $MPI_BIN_PATH
        done
	MPI_BIN_PATH=`pwd`
	cd $mydir	
	MPICC="$mympicc"

],
[C++], [
        AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_PATH_PROGS(MPICXX, 
		mpic++ mpiCC mpCC hcp mpxlC mpxlC_r cmpic++ sxmpic++, nocommand)
	if test x"$MPI_BIN_PATH" = "xnocommand"; then
     	     AC_MSG_ERROR([Can't find MPICXX])	
	fi		
	CXX="$MPICXX"
        AC_PROG_CXX()
	mydir=`pwd`
	mympicxx="$MPICXX"
	MPI_BIN_PATH=`dirname "$MPICXX"`
	cd $MPI_BIN_PATH
        while test -h "$MPICXX"; do
             MPICXX=`readlink "$MPICXX"`
             MPI_BIN_PATH=`dirname "$MPICXX"`
	     cd $MPI_BIN_PATH
        done
	MPI_BIN_PATH=`pwd`	
	cd $mydir	
	MPICXX="$mympicxx"
],
[Fortran 77], [
        AC_ARG_VAR(MPIF77,[MPI Fortran compiler command])
	AC_PATH_PROGS(MPIF77,
		mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c sxmpif77,
	     	nocommand)
	if test x"$MPI_BIN_PATH" = "xnocommand"; then
     	     AC_MSG_ERROR([Can't find MPIF77])	
	fi		
	F77="$MPIF77"
        AC_PROG_F77()
	mydir=`pwd`
	mympif77="$MPIF77"
	MPI_BIN_PATH=`dirname "$MPIF77"`
	cd $MPI_BIN_PATH
        while test -h "$MPIF77"; do
             MPIF77=`readlink "$MPIF77"`
             MPI_BIN_PATH=`dirname "$MPIF77"`
	     cd $MPI_BIN_PATH
        done
	MPI_BIN_PATH=`pwd`	
	cd $mydir	
	MPIF77="$mympif77"
],
[Fortran 90], [
        AC_ARG_VAR(MPIF90,[MPI Fortran compiler command])
        AC_PATH_PROGS(MPIF90,
		mpif90 hf90 mpxlf mpf90 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c sxmpif90,
	     	nocommand)
	if test x"$MPI_BIN_PATH" = "xnocommand"; then
     	     AC_MSG_ERROR([Can't find MPIF90])	
	fi		
	F90="$MPIF90"
        AC_PROG_F90()
	mydir=`pwd`
	mympif90="$MPIF90"
	MPI_BIN_PATH=`dirname "$MPIF90"`
	cd $MPI_BIN_PATH
        while test -h "$MPIF90"; do
             MPIF90=`readlink "$MPIF90"`
             MPI_BIN_PATH=`dirname "$MPIF90"`
	     cd $MPI_BIN_PATH
        done
	MPI_BIN_PATH=`pwd`	
	cd $mydir	
	MPIF90="$mympif90"
])

AC_SUBST(MPI_BIN_PATH)

dnl get mpi root path
mpi_bin_base=`basename "$MPI_BIN_PATH"`
if test x"$mpi_bin_base" = "xbin"; then
   MPI_ROOT_PATH=`AS_DIRNAME(["$MPI_BIN_PATH"])`
else
   AC_MSG_ERROR([Can't determine MPI_ROOT_PATH])	
fi	
AC_SUBST(MPI_ROOT_PATH)

MPI_LIB_PATH="$MPI_ROOT_PATH/lib"
dnl if ! test -d "$MPI_LIB_PATH" ; then
dnl    AC_MSG_ERROR([Can't determine MPI_LIB_PATH])
dnl fi


AC_LANG_CASE([C], [
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
],
[C++], [
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
])


])dnl AC_MPI
dnl

AC_DEFUN([AC_MPI_VERSION],[

echo -n "checking MPI version... "
MPI_VERSION=""
for file in `ls SUPPORTED_MPI_VERSIONS/*`; do
    if test x"$MPI_VERSION" = "x"; then
        . $file	
  	MPI_VERSION=`check_version`
    fi
done

dnl HP MPI
mpirun="$MPI_BIN_PATH/mpirun"
resh=`"$mpirun" -version 2>&1 | grep 'HP MPI'`
if test $? -eq 0; then
  ver=`echo $resh | sed 's/.*HP MPI \(.*\)/\1/' | sed 's/ /_/g'`
  if test $? -eq 1; then 
    ver="UNKNOWN"
  fi
  device="UNKNOWN"
  MPI_VERSION="HP_MPI-"$device"-"$ver
fi
dnl end HP MPI
dnl IBM MPI 
dnl POE
dnl Quadrics
dnl SUN
dnl

echo $MPI_VERSION
AC_SUBST(MPI_VERSION)

])

AC_DEFUN([AC_PROG_PERL_VERSION],[dnl
# Make sure we have perl
if test -z "$PERL"; then
   AC_CHECK_PROG(PERL,perl,perl)
fi

# Check if version of Perl is sufficient
ac_perl_version="$1"

if test x"$PERL" != "x"; then
  AC_MSG_CHECKING(for perl version greater than or equal to $ac_perl_version)
  # NB: It would be nice to log the error if there is one, but we cannot rely
  # on autoconf internals
  $PERL -e "use $ac_perl_version;" > /dev/null 2>&1
  if test $? -ne 0; then
    AC_MSG_RESULT(no);
    $3
  else
    AC_MSG_RESULT(ok);
    $2
  fi
else
  AC_MSG_WARN(Can't find PERL)
fi
])dnl


AC_DEFUN([AC_FS],[
CONFIGDIR=`pwd`
AC_SUBST(CONFIGDIR)
AC_ARG_VAR(SCRDIR,[Working directory])
if test x"$SCRDIR" = "x"; then
   SCRDIR=`pwd`
fi
])


AC_DEFUN([AC_PA],[
AC_ARG_VAR(PROC_ARRAY,[Number of used processors, e.g "2 4 8 16"])
if test x"$PROC_ARRAY" = "x"; then
   PROC_ARRAY="2"
fi
])



AC_DEFUN([AC_BATCH],[

AC_ARG_VAR(BATCH_SYSTEM,[See subdirectory SUPPORTED_BATCH_SYSTEMS for supported versions])
AC_ARG_VAR(QSUB,[qsub command])

dnl batch system

echo "checking batch system... "

AC_PATH_PROGS(QSUB, qsub bsub job_submit srun ll_submit condor_submit, nocommand)
if test x"$QSUB" != "xnocommand"; then
    if test x"$BATCH_SYSTEM" = "x"; then	
    	BATCH_BIN_PATH=`dirname "$QSUB"`
    	myqsub="$QSUB"
    	mydir=`pwd`
	cd $BATCH_BIN_PATH
    	while test -h "$QSUB"; do
        	QSUB=`readlink "$QSUB"`
        	BATCH_BIN_PATH=`dirname "$QSUB"`
        	cd $BATCH_BIN_PATH
    	done
    	BATCH_BIN_PATH=`pwd`	
    	cd $mydir	
    	QSUB="$myqsub"
    	QSUB_BASENAME=`basename "$QSUB"`    
    	if test x"$BATCH_SYSTEM" = "x"; then	
        	for file in `ls SUPPORTED_BATCH_SYSTEMS/*`; do
            	if test x"$BATCH_SYSTEM" = "x"; then
                	. $file	
                	BATCH_SYSTEM=`check_version`
            	fi
        	done
    	fi
    	if test x"$BATCH_SYSTEM" = "x"; then	
        	AC_MSG_ERROR([Can't determine batch system])	
    	fi
    fi
    file=SUPPORTED_BATCH_SYSTEMS/"$BATCH_SYSTEM"
    if test -f "$file"; then
        . $file
	set_version
        if test x"$QSUB" = "x"; then			
	    QSUB="$QSUB_BASENAME"
        fi
    else
        AC_MSG_ERROR([Can't determine batch system])	
    fi
else
    if test x"$BATCH_SYSTEM" = "x"; then
        BATCH_SYSTEM="NONE"
	QSUB="/bin/sh"
    fi
    file=SUPPORTED_BATCH_SYSTEMS/"$BATCH_SYSTEM"
    if test -f "$file"; then
        . $file
        set_version
        if test x"$QSUB" = "x"; then			
	    QSUB="$QSUB_BASENAME"
        fi
    else
        AC_MSG_ERROR([Can't determine batch system])	
    fi
fi

echo "checking batch system... "$BATCH_SYSTEM

AC_SUBST(QSUB)
AC_SUBST(QCONF)
AC_SUBST(BL1)
AC_SUBST(BL2)
AC_SUBST(BL3)
AC_SUBST(BL4)
AC_SUBST(BL5)
AC_SUBST(BL6)



dnl
])


AC_DEFUN([AC_NODE_NUM_PROCS],[
AC_ARG_VAR(NODE_NUM_PROCS,[MPI processors per node])
if test x"$NODE_NUM_PROCS" = "x"; then
  NODE_NUM_PROCS=1
fi
])


AC_DEFUN([AC_QUEUE],[
AC_ARG_VAR(QUEUE,[Batch queue name])
if test x"$QUEUE" = "x"; then
 QUEUE=""
fi
])

















