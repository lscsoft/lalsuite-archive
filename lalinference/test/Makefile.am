AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinference.la

if LALSUPPORT
if LALINSPIRAL
if LALXML

XMLPRG = \
	LALInferenceXMLTest

endif

LALInferenceTest_SOURCES = \
	LALInferenceTest.c \
	LALInferenceTest.h

LALInferencePriorTest_SOURCES = \
	LALInferencePriorTest.c
	LALInferenceTest.h

LALInferenceLikelihoodTest_SOURCES = \
	LALInferenceLikelihoodTest.c
	LALInferenceTest.h

LALInferenceProposalTest_SOURCES = \
	LALInferenceProposalTest.c
	LALInferenceTest.h

LALInferenceInterpsTest_SOURCES = \
	LALInferenceInterps.c \
	LALInferenceInterps.h

check_PROGRAMS = \
	LALInferenceTest \
	LALInferencePriorTest \
	LALInferenceLikelihoodTest \
	LALInferenceProposalTest \
	LALInferenceInterpsTest
	$(XMLPRG)

TESTS = \
	LALInferenceTest \
	LALInferencePriorTest \
	LALInferenceLikelihoodTest \
	LALInferenceInterpsTest \
	$(XMLPRG)

endif
endif

EXTRA_DIST =

CLEANFILES = \
	test_vot.xml \
	*.out \
	*.dat \
	*.xml

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
