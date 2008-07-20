MAINTAINERCLEANFILES = Makefile.in \
	aclocal.m4 config.guess config.h.in config.sub configure \
	depcomp install-sh ltmain.sh missing

EXTRA_DIST = autogen.sh

SUBDIRS = src doc

dist_doc_DATA = README COPYING AUTHORS ChangeLog INSTALL NEWS

if COND_PKGCONFIG
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = blitzwave.pc
endif
