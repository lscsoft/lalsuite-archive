PACKAGE=\"powerflux\"
VERSION=\"1.1.20\"
LALROOT=./../../../../
LALINCLUDE=-I$(LALROOT)/include/
LALLDFLAGS=-L$(LALROOT)/lib -L$(LALROOT)/support-lib -llalsupport -llal -lm
LDFLAGS=$(LALLDFLAGS) -lpng -lz
#OPTFLAGS=-O2 -funroll-loops 
OPTFLAGS=-O6 -mmmx -msse -mcpu=i686
OPTFLAGS=-O6 -mcpu=i686 -funroll-loops
#OPTFLAGS=-O6 -msse2 -mcpu=i686 -funroll-loops
CFLAGS=-Wall -g $(OPTFLAGS) -DPACKAGE=$(PACKAGE) -DVERSION=$(VERSION) \
	$(LALINCLUDE)
CC=gcc $(CFLAGS)

OBJS=cmdline.o rastermagic.o hookup.o fine_grid.o lines.o grid.o powerflux.o  \
	png.o polarization.o intervals.o statistics.o

# 199.75 hz
BAND=359650
# 200 hz
#BAND=360000
# 352 hz
#BAND=633600
# 650 Hz
#BAND=1170000
# 650.25 hz
#BAND=1170450

#
# Do not build anything by default
#
none:

all: powerflux


test1: powerflux
	time ./powerflux -i /nfs/node2-storage/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO

test1a: powerflux
	time ./powerflux -i /nfs/node2-storage-2/LIGO/S2/S2.H1.op/power.H1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO

test1b: powerflux
	time ./powerflux -i /nfs/node2-storage/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/nfs/node2-storage/LIGO/detresponse \
		--detector=LLO

test1c: powerflux
	time ./powerflux -i /nfs/node2-storage-2/LIGO/S2/S2.H1.op/power.H1. -f $(BAND) -n 501 \
		--ephemeris-path=/nfs/node2-storage/LIGO/detresponse \
		--detector=LHO

test2: test_rm
	./test_rm

test3: powerflux
	time ./powerflux --input-format=Power -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0 --fine-factor=5 \
		--subtract-background=1 --ks-test=1
#		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.1 \
		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \

test4: powerflux
	valgrind --tool=cachegrind ./powerflux -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0

test4b: powerflux.hydra
	valgrind --tool=memcheck --error-limit=no ./powerflux.hydra --input-format=Power -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0 \
		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.1 \
		 --only-large-cos=0.34
#		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \

test5: powerflux
	time ./powerflux -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0 --filter-lines=1 \
		--fake-strain=0.3e-22 --config=test_conf --resolution-ratio=1.0 \
		--fake-orientation=0.5 --fake-ra=5.785685 --fake-dec=-1.563849

test.H1.inj.1: powerflux
	time ./powerflux -i /home/volodya/LIGO/S2/inj/S2.H1.op.inj/power.H1. -f 2302000 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --side-cut=300 --nbins=501 --do-cutoff=0 --filter-lines=0 \
		--patch-N-RA=120 --patch-N-DEC=90

test.L1.inj.1: powerflux
	time ./powerflux -i /home/volodya/LIGO/S2/inj/S2.L1.op.inj/power.L1. -f 2302000 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --side-cut=300 --nbins=501 --do-cutoff=0 --filter-lines=0 \
		--patch-N-RA=120 --patch-N-DEC=90
		


S3.H1.silver: powerflux
	time ./powerflux -i /home/volodya/LIGO/S3/S3.H1.op/power.H1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO
	

WL=-Wl,--rpath -Wl,/home/volodya/LIGO/LAL/support-lib/.libs -Wl,--rpath -Wl,/home/volodya/LIGO/LAL/lib/.libs -Wl,--rpath -Wl,/usr/local/lib -Wl,--rpath -Wl,/usr/local/lib

powerflux: $(OBJS)
	$(CC) $(OBJS) -o powerflux -static $(LDFLAGS) -lgsl -lgslcblas
	ls -l powerflux

powerflux.hydra: $(OBJS)
	(export PATH=/opt/condor/bin:$$PATH ; \
	 export CONDOR_CONFIG=/opt/condor-6.6.5/etc/condor_config ; \
	 condor_compile $(CC) $(OBJS) -o powerflux.hydra -static $(LDFLAGS) -lgsl -lgslcblas \
	 )
	ls -l powerflux.hydra
	

test_rm: test_rm.c rastermagic.o
	$(CC) test_rm.c rastermagic.o -lm -o test_rm

cmdline.c: powerflux.ggo
	gengetopt  --conf-parser < powerflux.ggo

clean:
	rm -f core *.bck *.o

upload: powerflux.hydra
	scp -S globus-ssh -C powerflux.hydra hydra.phys.uwm.edu:PowerFlux/powerflux
	
upload.g: powerflux.hydra
	scp -C powerflux.hydra gallatin.physics.lsa.umich.edu:PowerFlux/powerflux
	
