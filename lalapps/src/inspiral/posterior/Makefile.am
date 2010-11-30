SUBDIRS = SPINspiral

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALINSPIRAL
if LALFRAME
if LALMETAIO
metaio_progs = lalapps_inspnest lalinference_test
lalapps_inspnest_SOURCES = lalnest.c LALInspiralMCMC.c LALInspiralMCMCUser.c nest_calc.c priors.c MultiNest_calc.c nest_calc.h LALInspiralMCMC.h LALInspiralMCMCUser.h priors.h MultiNest_calc.h
lalapps_inspnest_LDADD = -L$(top_builddir)/src/inspiral/posterior/MultiNest_v2.8 -lnest3 -llapack -lblas -lgfortranbegin -lgfortran -lm
endif
lalinference_test_SOURCES = LALInference.c LALInferenceReadData.c LALInferenceTemplate.c InferenceTest.c
frame_progs = lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = followupMcmc.c
endif
endif

bin_PROGRAMS = $(metaio_progs) $(frame_progs)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps

if HAVE_PYTHON
bin_SCRIPTS = lalapps_combine_evidence lalapps_coherence_test lalapps_merge_nested_sampling_runs lalapps_nest_pipe
pkgpython_DATA = nest_utils.py combine_evidence.py
endif

EXTRA_DIST = \
	combine_evidence.in \
	coherence_test.in \
	merge_nested_sampling_runs.in \
	nest_pipe.in \
	nest_utils.py \
	combine_evidence.py

python_config.sed: Makefile
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_combine_evidence: $(srcdir)/combine_evidence.in python_config.sed
	sed -f python_config.sed $(srcdir)/combine_evidence.in > lalapps_combine_evidence
	chmod a+x lalapps_combine_evidence

lalapps_coherence_test: $(srcdir)/coherence_test.in python_config.sed
	sed -f python_config.sed $(srcdir)/coherence_test.in > lalapps_coherence_test
	chmod a+x lalapps_coherence_test

lalapps_merge_nested_sampling_runs: $(srcdir)/merge_nested_sampling_runs.in python_config.sed
	sed -f python_config.sed $(srcdir)/merge_nested_sampling_runs.in >lalapps_merge_nested_sampling_runs
	chmod a+x lalapps_merge_nested_sampling_runs

lalapps_nest_pipe: $(srcdir)/nest_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/nest_pipe.in >lalapps_nest_pipe
	chmod a+x lalapps_nest_pipe

CLEANFILES = $(bin_SCRIPTS) python_config.sed

BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
