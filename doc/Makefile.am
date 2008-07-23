MAINTAINERCLEANFILES = Makefile.in api/*

EXTRA_DIST = index.xhtml packagedoc.css api/index.html api examples

# apidir = $(docdir)/api
# api_DATA = api/index.html api/*

# examplesdir = $(docdir)/examples
# examples_DATA = examples/*

api/index.html: $(top_srcdir)/src/*.cxx $(top_srcdir)/src/*.h
	doxygen
