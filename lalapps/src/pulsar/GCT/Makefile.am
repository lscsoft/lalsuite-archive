MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic -I$(top_srcdir)/src/pulsar/EinsteinAtHome -I$(top_srcdir)/src/pulsar/HoughFstat

bin_PROGRAMS = lalapps_HierarchSearchGCT
check_PROGRAMS =
EXTRA_PROGRAMS = eah_HierarchSearchGCT eah_HierarchSearchGCT_SSE2 eah_HierarchSearchGCT_SSE2_NOOPT eah_HierarchSearchGCT_SSE2_INTRIN \
	         lalapps_HierarchSearchGCT_SSE2 lalapps_HierarchSearchGCT_SSE2_NONC

EAHSources = hs_boinc_extras.c hs_boinc_options.cpp win_lib.cpp
HSGCTSources = HierarchSearchGCT.c HierarchSearchGCT.h \
                          GCTtoplist.c GCTtoplist.h \
                          HeapToplist.c \
                          RecalcToplistStats.c RecalcToplistStats.h

lalapps_HierarchSearchGCT_SOURCES = $(HSGCTSources)

lalapps_HierarchSearchGCT_SSE2_SOURCES = $(HSGCTSources)
lalapps_HierarchSearchGCT_SSE2_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalapps_HierarchSearchGCT_SSE2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

lalapps_HierarchSearchGCT_SSE2_NONC_SOURCES = $(HSGCTSources)
lalapps_HierarchSearchGCT_SSE2_NONC_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT -DEXP_NO_NUM_COUNT
lalapps_HierarchSearchGCT_SSE2_NONC_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse


eah_HierarchSearchGCT_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_SSE2_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_SSE2_NOOPT_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_SSE2_INTRIN_SOURCES = $(HSGCTSources) $(EAHSources)

eah_HierarchSearchGCT_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT
eah_HierarchSearchGCT_SSE2_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT 
eah_HierarchSearchGCT_SSE2_NOOPT_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT -DEXP_NO_ASM
eah_HierarchSearchGCT_SSE2_INTRIN_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT -DEXP_USE_INTRIN

# this shouldn't harm - it just prefers libraries in the build directory over e.g. system ones
eah_HierarchSearchGCT_LDFLAGS = -L. $(AM_LDFLAGS)
eah_HierarchSearchGCT_SSE2_LDFLAGS = -L. $(AM_LDFLAGS)
eah_HierarchSearchGCT_SSE2_NOOPT_LDFLAGS = -L. $(AM_LDFLAGS)
eah_HierarchSearchGCT_SSE2_INTRIN_LDFLAGS = -L. $(AM_LDFLAGS)

# RELEASE_DEPS contains additional dependencies of the E@H build
# items in there should have rules listed below to create them
# the idea is to increase compatibility of E@H applications by linking these libs statically
eah_HierarchSearchGCT_DEPENDENCIES = $(RELEASE_DEPS)
eah_HierarchSearchGCT_SSE2_DEPENDENCIES = $(RELEASE_DEPS)
eah_HierarchSearchGCT_SSE2_NOOPT_DEPENDENCIES = $(RELEASE_DEPS)
eah_HierarchSearchGCT_SSE2_INTRIN_DEPENDENCIES = $(RELEASE_DEPS)

eah_HierarchSearchGCT_LDADD = $(LDADD) $(RELEASE_LDADD)
eah_HierarchSearchGCT_SSE2_LDADD = $(LDADD) $(RELEASE_LDADD)
eah_HierarchSearchGCT_SSE2_NOOPT_LDADD = $(LDADD) $(RELEASE_LDADD)
eah_HierarchSearchGCT_SSE2_INTRIN_LDADD = $(LDADD) $(RELEASE_LDADD)

libstdc++.a:
	cp -f `$(CC) -print-file-name=$@` $@

exchndl.o: exchndl.c
	$(CC) $(AM_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

erp_execinfo_plus.o: erp_execinfo_plus.c
	$(CC) $(AM_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


# how to compile CUDA files (.cu)
NVCC ?= nvcc
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@

EXTRA_DIST = testGCT.sh testGCT_inject.sh

TESTS = testGCT.sh testGCT_inject.sh

CLEANFILES =
