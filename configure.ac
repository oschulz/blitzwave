# Initialize:

AC_INIT([blitzwave], [0.7.1])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Checks for programs:

AC_PROG_CXX
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

AC_CHECK_PROGS(DOXYGEN, doxygen, false)
AM_CONDITIONAL([COND_DOXYGEN], [test "$DOXYGEN" != "false"])

AC_CHECK_PROGS(PKGCONFIG, pkg-config, false)
AM_CONDITIONAL([COND_PKGCONFIG], [test "$PKGCONFIG" != "false"])

# Checks for libraries:

REQ_PKGS=""
REQ_LIBS="blitz"

for LIBNAME in $REQ_LIBS; do
	LIBCONFIG="${LIBNAME}-config"
	echo "SEARCHING FOR \"${LIBCONFIG}\" ..."
	AC_MSG_CHECKING(${PKGCONFIG} ${LIBNAME})
	if ${PKGCONFIG} ${LIBNAME}; then
		AC_MSG_RESULT(yes)
		REQ_PKGS="$REQ_PKGS ${LIBNAME}"
		LIBS="`${PKGCONFIG} ${LIBNAME} --libs` $LIBS"
	else
		AC_MSG_RESULT(no)
		AC_MSG_CHECKING(${LIBCONFIG})
		if which "${LIBCONFIG}" > /dev/null 2> /dev/null ; then
			AC_MSG_RESULT(yes)
			LIB_CFLAGS="$LIB_CFLAGS `${LIBCONFIG} --cflags`"
			NEW_LIBS="`${LIBCONFIG} --libs`"
			LIB_LIBS="${NEW_LIBS} $LIB_LIBS"
			LIBS="${NEW_LIBS} $LIBS"
		else
			AC_MSG_RESULT(no)
			NEW_LIBS=""
			AC_CHECK_LIB(${LIBNAME}, main, NEW_LIBS="-l${LIBNAME}", [AC_MSG_ERROR([Could not find ]${LIBNAME}[ library!])])
			LIB_LIBS="${NEW_LIBS} $LIB_LIBS"
			LIBS="${NEW_LIBS} $LIBS"
		fi
	fi
done

AC_SUBST(LIB_CFLAGS)
AC_SUBST(LIB_LIBS)
AC_SUBST(REQ_PKGS)


# Output:

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
	Makefile
	blitzwave.pc
	src/Makefile
	doc/Makefile doc/Doxyfile
])

AC_OUTPUT
