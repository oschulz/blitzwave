dnl -*- mode: autoconf -*- 
dnl
dnl Autoconf macro resolve C/C++ library dependencies
dnl
dnl  CXX_PROJECT_DEPS([DEPENDENCIES])
dnl
dnl The macro defines the following substitution variables
dnl
dnl    DEP_CFLAGS         CFLAGS for depencencies, also added to CFLAGS
dnl    DEP_LIBS           LIBS for dependencies, also added to LIBS


AC_DEFUN([CXX_PROJECT_DEPS],
[
	AC_REQUIRE([CXX_PROJECT_INIT])

	DEPS="$1"

	PKGCONFIG_DEPS=""
	for DEPNAME in $DEPS; do
		echo "SEARCHING FOR \"${DEPNAME}\" ..."
		AC_MSG_CHECKING(${PKGCONFIG} ${DEPNAME})
		if ${PKGCONFIG} ${DEPNAME}; then
			AC_MSG_RESULT(yes)
			PKGCONFIG_DEPS="$PKGCONFIG_DEPS ${DEPNAME}"
		else
			AC_MSG_RESULT(no)
			LIBCONFIG="${DEPNAME}-config"
			AC_MSG_CHECKING(${LIBCONFIG})
			if which "${LIBCONFIG}" > /dev/null 2> /dev/null ; then
				AC_MSG_RESULT(yes)
				DEP_CFLAGS="${DEP_CFLAGS} "`"${LIBCONFIG}" --cflags`
				DEP_LIBS="`${LIBCONFIG} --libs` ${DEP_LIBS}"
			else
				AC_MSG_RESULT(no)
				NEW_LIBS=""
				AC_CHECK_LIB(${DEPNAME}, main, NEW_LIBS="-l${DEPNAME}", [AC_MSG_ERROR([Could not find ]${DEPNAME}[ library!])])
				DEP_LIBS="${NEW_LIBS} ${DEP_LIBS}"
			fi
		fi
	done

	CFLAGS="$CFLAGS $DEP_CFLAGS"
	CXXFLAGS="$CXXFLAGS $DEP_CFLAGS"
	LIBS="$DEP_LIBS $LIBS"

	AC_SUBST(DEP_CFLAGS)
	AC_SUBST(DEP_LIBS)

	if test "${PKGCONFIG_DEPS}" != "" ; then
		for PKGNAME in ${PKGCONFIG_DEPS}; do
			if ! "${PKGCONFIG}" "${PKGNAME}"; then
				AC_MSG_ERROR([Package requirements not met - no package ']${PKGNAME}[' found])
			fi
		done
		CFLAGS="$CFLAGS `${PKGCONFIG} ${PKGCONFIG_DEPS} --cflags`"
		CXXFLAGS="$CXXFLAGS `${PKGCONFIG} ${PKGCONFIG_DEPS} --cflags`"
		LIBS="`${PKGCONFIG} ${PKGCONFIG_DEPS} --libs` $LIBS"
	fi

	AC_SUBST(PKGCONFIG_DEPS)

	SHLIBEXT="${shrext_cmds}"
	AC_SUBST(SHLIBEXT)
])


#
# EOF
#
