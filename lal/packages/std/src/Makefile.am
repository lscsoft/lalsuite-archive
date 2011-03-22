INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

if CUDA
MALLOCSRC = LALMallocCUDA.c
else
MALLOCSRC = LALMalloc.c
endif

noinst_LTLIBRARIES = libstd.la
libstd_la_SOURCES = \
	Grid.c \
	LALComplex.c \
	LALError.c \
	LALGSL.c \
	LALStd.c \
	LALString.c \
	LALVersion.c \
	StringConvert.c \
	StringToken.c \
	XLALError.c \
	$(MALLOCSRC)

EXTRA_DIST = \
	Grid.m4 \
	LALCreateDestroyGrid.m4 \
	LALPrimer.c

EXTRA_TEXSOURCES = LALPrimer.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = \
	Grid.c

dvi-local: .dvi-dep

Grid.c: Grid.m4 LALCreateDestroyGrid.m4
	cd $(srcdir) && @M4@ Grid.m4 > Grid.c

@DVIMKFILE@


# cuda build rules
#include ${top_srcdir}/misc/cuda.common

if CUDA
SUFFIXES = .cu .cpp
.cu.cpp:
	@nvcc --cuda $(INCLUDES) $(CPPFLAGS) $(NVCC_CFLAGS) --output-file $@ $<
endif

