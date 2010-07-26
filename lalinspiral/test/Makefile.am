AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

if LALSUPPORT
check_PROGRAMS = \
	BCVTemplates \
	BCVTemplatesFlatMesh \
	BasicInjectTest \
	CoarseTest \
	CoarseTest2 \
	FindChirpBCVSpinTest \
	FindChirpBankVetoTest \
	FindChirpTDTest \
	FindChirpACTDTest \
        FisherTest \
	GenerateInspiralWaveform \
        GeneratePPNAmpCorConsistencyTest \
	GeneratePPNAmpCorInspiralTest \
	GeneratePPNInspiralTest \
	GetOrientationEllipse \
	InjectionInterfaceTest \
	InspiralBCVSpinBankTest \
	InspiralSpinBankTest \
	LALInspiralSpinningBHBinariesTest \
	LALInspiralTest \
	LALSTPNWaveformTest \
	MetricTest \
	MetricTestBCV \
	MetricTestPTF \
	PNTemplates \
	RandomInspiralSignalTest
TESTS = $(check_PROGRAMS)
FindChirpTDTest_SOURCES = FindChirpTDTest.c FindChirpTDTest.h
FindChirpACTDTest_SOURCES = FindChirpACTDTest.c FindChirpACTD.h FindChirpTDTest.h
endif

# non-building test
# BCVSpinTemplates
# ChirpSpace
# FilterTest
# FindChirpBCVTest
# FindChirpSPTest
# FindChirpTmpltTest
# InspiralSpinBankwNDTemplateBankTest
# NRWaveInjectTest
# RingTest
# TestRingdown

# failing test
# LALHybridTest
# NRWaveIOTest
# SpaceCovering

# requires lalstochastic
# CoherentInspiralFilterTest
# CoherentInspiralInputTest

EXTRA_DIST = \
	NRWaveIOTest.cfg \
	NRWaveIOTest.data \
	RandomInspiralSignalTest.par \
	injection.xml

CLEANFILES = \
	*.out \
	*.dat

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
