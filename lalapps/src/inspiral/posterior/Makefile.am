SUBDIRS = SPINspiral

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALINSPIRAL
if LALFRAME

if LALMETAIO
metaio_progs = lalapps_inspnest lalinference_test
lalapps_inspnest_SOURCES = lalnest.c LALInspiralMCMC.c LALInspiralMCMCUser.c nest_calc.c priors.c MultiNest_calc.c nest_calc.h LALInspiralMCMC.h LALInspiralMCMCUser.h priors.h MultiNest_calc.h
lalapps_inspnest_LDADD = -llapack -L$(top_builddir)/src/inspiral/posterior/MultiNest_v2.8 -lnest3
endif
lalinference_test_SOURCES = LALInference.c LALInferenceReadData.c LALInferenceTemplate.c InferenceTest.c

frame_progs = lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = followupMcmc.c
endif
endif

bin_PROGRAMS = $(metaio_progs) $(frame_progs)
