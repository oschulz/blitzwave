MAINTAINERCLEANFILES = Makefile.in

EXTRA_DIST = testimage.pnm

lib_LTLIBRARIES = libblitzwave.la

libblitzwave_la_SOURCES = arrayTools.cxx Wavelet.cxx WaveletDecomp.cxx
libblitzwave_la_headers = arrayTools.h WaveletDecomp.h Wavelet.h

include_HEADERS = $(libblitzwave_la_headers)

noinst_PROGRAMS = arrayTools-test WaveletDecomp-test Wavelet-test

arrayTools_test_SOURCES = arrayTools-test.cxx
arrayTools_test_LDADD = libblitzwave.la
arrayTools_test_LDFLAGS = -static

WaveletDecomp_test_SOURCES = WaveletDecomp-test.cxx
WaveletDecomp_test_LDADD = libblitzwave.la
WaveletDecomp_test_LDFLAGS = -static

Wavelet_test_SOURCES = Wavelet-test.cxx
Wavelet_test_LDADD = libblitzwave.la
Wavelet_test_LDFLAGS = -static
