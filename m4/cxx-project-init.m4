dnl -*- mode: autoconf -*- 
dnl
dnl Autoconf macro initialize a basic C++ project with a few custom options
dnl Synopsis:
dnl
dnl  CXX_PROJECT_INIT
dnl


AC_DEFUN([CXX_PROJECT_INIT],
[
	# Workaround for older autoconf versions (< 2.60)
	if test "x$docdir" == "x"; then
		docdir='${datadir}/doc/${PACKAGE_TARNAME}'
		AC_SUBST(docdir)
	fi

	# Checks for tools:

	AC_CHECK_PROGS(GREP, grep, false)
	AC_CHECK_PROGS(SED, sed, false)

	AC_CHECK_PROGS(DOXYGEN, doxygen, false)
	AM_CONDITIONAL([COND_DOXYGEN], [test "$DOXYGEN" != "false"])

	AC_CHECK_PROGS(PKGCONFIG, pkg-config, false)
	AM_CONDITIONAL([COND_PKGCONFIG], [test "$PKGCONFIG" != "false"])
])


#
# EOF
#
