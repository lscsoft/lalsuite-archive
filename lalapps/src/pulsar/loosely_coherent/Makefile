PACKAGE=\"loosely_coherently_coherent\"
VERSION=0.20.0-64
LALROOT=/opt/lalsuite/
LALINCLUDE=-I$(LALROOT)/include/
LALLDFLAGS=-L$(LALROOT)/lib -L$(LALROOT)/support-lib -L/opt/zlib/lib -L/opt/png/lib -llalsupport -llalpulsar -llal -lm -lz
LDFLAGS=$(LALLDFLAGS) -L. -lpng -lfftw3 -lfftw3f
#OPTFLAGS=-O2 -funroll-loops 
OPTFLAGS=-O6 -mmmx -msse -mcpu=i686
OPTFLAGS=-O6 -march=i686 -funroll-loops
OPTFLAGS=-O6 -msse2 -march=i686 -funroll-loops
OPTFLAGS=-O6 -march=opteron -funroll-loops -fno-stack-protector
OPTFLAGS=-O -march=opteron -fno-stack-protector
OPTFLAGS=-O6 -msse2 -funroll-loops -fno-stack-protector
#OPTFLAGS=-O6 -mtune=opteron -mfpmath=sse -funroll-loops -fno-stack-protector -U_FORTIFY_SOURCE
OPTFLAGS=-O6 -mtune=core2 -funroll-loops -fno-stack-protector -msse3 -U_FORTIFY_SOURCE
CFLAGS=-Wall -g $(OPTFLAGS) -DPACKAGE=$(PACKAGE) -DVERSION=\"$(VERSION)\" \
	$(LALINCLUDE) -fPIC
CC=gcc $(CFLAGS)

COMMON_OBJS=cmdline.o rastermagic.o lines.o polarization.o \
	png.o intervals.o dataset.o skymarks.o hookup.o statistics.o grid.o bessel.o context.o fft_stats.o

OBJS=$(COMMON_OBJS) loosely_coherent1.o

FILES=no_marks.txt sky_marks.txt all_sky_marks_*.txt

# 50 hz
#BAND=90000
# 140.5 hz
#BAND=252900
# 149 hz
#BAND=268200
# 199.75 hz
#BAND=359650
# 200 hz
#BAND=360000
# 400 Hz
# 352 hz
#BAND=633600
# 650 Hz
#BAND=1170000
# 650.25 hz
#BAND=1170450

#DATASET=/home/volodya/LIGO/S4/Injections/sfts/108.dst
DATASET=$(HOME)/LIGO/SFT-DATA/S5/H1_395_415.dst
BAND=1600
#DATASET=random.dst
RA=0.4687356
DEC=-1.257051

RA=3.14
DEC=0

# DATASET=$(HOME)/LIGO/SFT-DATA/Gaussian/H1_400_410.dst
# BAND=1626
# #DATASET=random.dst
# RA=1.0
# DEC=0.1
# F0=406.527
# F1=-1.35e-10
# REFTIME=850000000


#
# Do not build anything by default
#
none:

all: loosely_coherent1



test3c: loosely_coherent1
	./loosely_coherent1 --dataset=$(DATASET) -f $(BAND) -n 1 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalsuite/lalpulsar/test/detresponse \
		--filter-lines=1 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lalsuite/lalpulsar/test//earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lalsuite/lalpulsar/test//sun05-09.dat  \
		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.0001 --focus-f0=400.01 --focus-f0-delta=0.001 \
		--dump-data=input.dat \
		--write-dat=NONE \
		--write-png=.*$$ 

test3d: loosely_coherent1
	./loosely_coherent1 --config=config.3d
	head 1626-406.500000/data.log
#
#		--fake-linear --fake-ra=3.141593 --fake-dec=0.006233 --fake-freq=199.85 --fake-strain=1e-20 \
#		--fake-spindown=1e-9\
#		 --band-axis=explicit\(0.5,0.3,0.1\)


test3e: loosely_coherent1
	./loosely_coherent1 --config=config.3e
	head 1626-406.500000/data.log

test3j: loosely_coherent1
	./loosely_coherent1 --config=config.3j
	head 1626-406.500000/data.log

test3jmt: loosely_coherent1-mt
	./loosely_coherent1-mt --config=config.3j
	head 1626-406.500000/data.log

test_make_dst: loosely_coherent1
	./loosely_coherent1 --config=config.make_dst
	head 1626-406.500000/data.log

test2cl3e: loosely_coherent2cl
	./loosely_coherent2cl --config=config.3e
	head 1626-406.500000/data.log

test2cl3j: loosely_coherent2cl
	./loosely_coherent2cl --config=config.3j
	head 1626-406.500000/data.log

dat2skymap: grid.o rastermagic.o dat2skymap.o png.o
	$(CC) grid.o rastermagic.o dat2skymap.o png.o -o dat2skymap -static $(LDFLAGS)



loosely_coherent1: $(OBJS) jobs.o util.o
	$(CC) $(OBJS) jobs.o util.o -o loosely_coherent1 -static $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l loosely_coherent1

loosely_coherent1-mt: $(OBJS) jobs-mt.o util.o
	$(CC) $(OBJS) jobs-mt.o util.o -o loosely_coherent1-mt -static $(LDFLAGS) -lgsl -lgslcblas -lm -lpthread -lz
	ls -l loosely_coherent1-mt

loosely_coherent1.dyn: $(OBJS) jobs.o util.o
	$(CC) $(OBJS) jobs.o util.o -o loosely_coherent1.dyn $(LDFLAGS) -lgsl -lgslcblas -lz
	ls -l loosely_coherent1.dyn

loosely_coherent1-condor: $(OBJS) jobs.o util-condor.o
	 condor_compile $(CC) $(OBJS) jobs.o util-condor.o -o loosely_coherent1-condor -static $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l loosely_coherent1-condor

loosely_coherent2cl: $(OBJS) jobs.o util.o ocl1.o loosely_coherent2cl.o
	$(CC) $(COMMON_OBJS) loosely_coherent2cl.o jobs.o util.o ocl1.o -lOpenCL -o loosely_coherent2cl  $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l loosely_coherent2cl

ocl1.o: ocl1.c ocl1.h
	$(CC) $(CFLAGS) -I/usr/local/cuda/include -c ocl1.c

jobs-mt.o: jobs.c jobs.h
	$(CC) $(CFLAGS) -DLINUX_THREADS -c jobs.c -o jobs-mt.o

util-condor.o: util.c util.h
	$(CC) $(CFLAGS) -DCONDOR -c util.c -o util-condor.o

test_rm: test_rm.c rastermagic.o
	$(CC) test_rm.c rastermagic.o -lm -o test_rm

cmdline.c: loosely_coherent.ggo
	gengetopt  --string-parser --conf-parser < loosely_coherent.ggo

clean:
	rm -f core *.bck *.o

upload: loosely_coherent1-condor loosely_coherent1-mt
	gsiscp -C loosely_coherent1-condor hydra.phys.uwm.edu:PowerFlux/loosely_coherent1-condor.$(VERSION)
	gsiscp -C loosely_coherent1-mt hydra.phys.uwm.edu:PowerFlux/loosely_coherent1-mt.$(VERSION)
	gsiscp -C $(FILES) hydra.phys.uwm.edu:PowerFlux/
	
upload.a: loosely_coherent1-condor loosely_coherent1-mt
	gsiscp -C loosely_coherent1-condor h2.atlas.aei.uni-hannover.de:PowerFlux/loosely_coherent1-condor.$(VERSION)
	gsiscp -C loosely_coherent1-mt h2.atlas.aei.uni-hannover.de:PowerFlux/loosely_coherent1-mt.$(VERSION)
	gsiscp -C $(FILES) h2.atlas.aei.uni-hannover.de:PowerFlux/
	
upload.g: loosely_coherent1-condor loosely_coherent1-mt
	scp -C loosely_coherent1-condor gallatin.physics.lsa.umich.edu:PowerFlux/loosely_coherent1-condor.$(VERSION)
	scp -C loosely_coherent1-mt gallatin.physics.lsa.umich.edu:PowerFlux/loosely_coherent1-mt.$(VERSION)
	scp -C $(FILES) gallatin.physics.lsa.umich.edu:PowerFlux/
	
upload.h: loosely_coherent1-condor loosely_coherent1-mt
	gsiscp -C loosely_coherent1-condor ldas-grid.ligo-wa.caltech.edu:PowerFlux/loosely_coherent1-condor.$(VERSION)
	gsiscp -C $(FILES) ldas-grid.ligo-wa.caltech.edu:PowerFlux/

upload.l: loosely_coherent1-condor loosely_coherent1-mt
	gsiscp -C loosely_coherent1-condor ldas-grid.ligo-la.caltech.edu:PowerFlux/loosely_coherent1-condor.$(VERSION)
	gsiscp -C $(FILES) ldas-grid.ligo-la.caltech.edu:PowerFlux/

upload.c: loosely_coherent1-condor loosely_coherent1-mt 
	gsiscp -C loosely_coherent1-condor ldas-grid.ligo.caltech.edu:PowerFlux/loosely_coherent1-condor.$(VERSION)
	gsiscp -C loosely_coherent1-mt ldas-grid.ligo.caltech.edu:PowerFlux/loosely_coherent1-mt.$(VERSION)
	gsiscp -C $(FILES) ldas-grid.ligo.caltech.edu:PowerFlux/

upload.m: loosely_coherent1-condor loosely_coherent1-mt
	gsiscp -C loosely_coherent1-condor vdergach@morgane.aei.mpg.de:PowerFlux/loosely_coherent1-condor.$(VERSION)
	gsiscp -C loosely_coherent1-mt vdergach@morgane.aei.mpg.de:PowerFlux/loosely_coherent1-mt.$(VERSION)
	gsiscp -C $(FILES) vdergach@morgane.aei.mpg.de:PowerFlux/
	
