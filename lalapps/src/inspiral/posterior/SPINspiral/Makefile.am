if LALFRAME
if LALMETAIO
if LALINSPIRAL
bin_PROGRAMS = lalapps_spinspiral
lalapps_spinspiral_SOURCES = \
	SPINspiral_main.c \
	SPINspiral_mcmc.c \
	SPINspiral_multinest.c \
	SPINspiral_templates.c \
	SPINspiral_routines.c \
	SPINspiral_data.c \
	SPINspiral_parameters.c \
	SPINspiral_signal.c \
	SPINspiral_lal.c \
	SPINspiral_3rdparty.c \
	SPINspiral.h remez.h
lalapps_spinspiral_LDADD = -L$(top_builddir)/src/inspiral/posterior/MultiNest_v2.8 -lnest3 -llapack -lblas -lgfortranbegin -lgfortran -lm
endif
endif
endif

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

clean-local :
	rm -f mcmc.output.??????.?? *-signal*.dat
