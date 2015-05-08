# Initialize:

AC_INIT([blitzwave], [0.8.0])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror])


# Libtool init and language selection

# Workaround for newer vs. older versions of Automake:
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AC_PROG_CXX

CXX_PROJECT_INIT


# Checks for dependencies

CXX_PROJECT_DEPS([blitz])



# Output

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
	Makefile
	blitzwave.pc
	blitzwave-config
	src/Makefile
	doc/Makefile doc/Doxyfile
])

AC_OUTPUT
