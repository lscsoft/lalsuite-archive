PACKAGE=\"powerflux\"
VERSION=1.4.5
LALROOT=./../../../../
LALROOT=/opt/lscsoft/lal/
LALINCLUDE=-I$(LALROOT)/include/
LALLDFLAGS=-L$(LALROOT)/lib -L$(LALROOT)/support-lib -llalsupport -llal -lm
LDFLAGS=$(LALLDFLAGS) -lpng -lz
#OPTFLAGS=-O2 -funroll-loops 
OPTFLAGS=-O6 -mmmx -msse -mcpu=i686
OPTFLAGS=-O6 -march=i686 -funroll-loops
OPTFLAGS=-O6 -msse2 -march=i686 -funroll-loops
OPTFLAGS=-O6 -msse2 -march=opteron -funroll-loops
OPTFLAGS=-O6 -msse2  -funroll-loops
CFLAGS=-Wall -g $(OPTFLAGS) -DPACKAGE=$(PACKAGE) -DVERSION=\"$(VERSION)\" \
	$(LALINCLUDE)
CC=gcc $(CFLAGS)

OBJS=cmdline.o rastermagic.o hookup.o fine_grid.o lines.o grid.o powerflux.o  \
	png.o polarization.o intervals.o statistics.o dataset.o candidates.o

# 50 hz
#BAND=90000
# 140.5 hz
BAND=252900
# 149 hz
#BAND=268200
# 199.75 hz
#BAND=359650
# 200 hz
#BAND=360000
# 352 hz
#BAND=633600
# 650 Hz
#BAND=1170000
# 650.25 hz
#BAND=1170450

#DATASET=/home/volodya/LIGO/S4/Injections/sfts/108.dst
DATASET=/home/volodya/LIGO/S4/Injections/sfts/58.dst
RA=0.4687356
DEC=-1.257051


#
# Do not build anything by default
#
none:

all: powerflux.hydra


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
	time ./powerflux --input-format=Power -i /home/volodya/LIGO/S3/S3.H1.op/power.H1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --three-bins=0 --do-cutoff=0 \
		--subtract-background=1 --ks-test=1 \
		--spindown-start=0.0 --spindown-step=1e-9 --spindown-count=2
#		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.1 \
		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \

test3a: powerflux
	time ./powerflux --input-format=Power -i /home/volodya/LIGO/S3/S3.H1.op/power.H1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --three-bins=0 --do-cutoff=0 \
		--subtract-background=0 --ks-test=1 --filter-lines=0 \
		--spindown-start=1e-9 \
		--fake-linear --fake-ra=3.141593 --fake-dec=0.006233 --fake-freq=199.85 --fake-strain=1e-20 \
		--fake-spindown=1e-9
#		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.1 \
		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \

test3b: powerflux
	./powerflux --input-format=Power --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --three-bins=0 --do-cutoff=0 \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0e-9 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat \

#		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.2 \

#		--fake-linear --fake-ra=3.141593 --fake-dec=0.006233 --fake-freq=199.85 --fake-strain=1e-20 \
#		--fake-spindown=1e-9
#		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \

test3c: powerflux
	./powerflux --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--three-bins=0 --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0 --dump-candidates=0 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.7 \
		--dump-data=input.dat \
		--write-dat=NONE \
		--write-png=ul\\.png$$ 
#
#		--fake-linear --fake-ra=3.141593 --fake-dec=0.006233 --fake-freq=199.85 --fake-strain=1e-20 \
#		--fake-spindown=1e-9\
#		 --band-axis=explicit\(0.5,0.3,0.1\)

test4: powerflux
	valgrind --tool=cachegrind ./powerflux -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0 --ks-test=1 --subtract-background=1

test4a: powerflux.dyn
	( export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LALROOT)/lib:$(LALROOT)/support-lib ; \
	valgrind --tool=memcheck  --error-limit=no --leak-check=full ./powerflux.dyn --input-format=Power -i /home/volodya/LIGO/S3/S3.H1.op/power.H1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --three-bins=0 --do-cutoff=0 \
		--subtract-background=1 --ks-test=1 \
		--spindown-start=0.0 --spindown-step=1e-9 --spindown-count=2 \
		--focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.03 \
		--write-dat=NONE --write-png=NONE \
	)

test4b: powerflux.hydra
	valgrind --tool=memcheck --error-limit=no ./powerflux.hydra --input-format=Power -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0\
		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.1 \
		 --only-large-cos=0.34
#		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \

test4c: powerflux.dyn
	( export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LALROOT)/lib:$(LALROOT)/support-lib ; \
	valgrind --tool=memcheck --error-limit=no  --leak-check=full ./powerflux.dyn --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--three-bins=0 --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0 --dump-candidates=0 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.7 \
		--write-dat=NONE \
		--write-png=ul\\.png\\$$ \
	)
#		 --focus-ra=2.2 --focus-dec=0.3 --focus-radius=0.1 \
		 --band-axis=explicit\(0.5,0.3,0.1\)
#		 --write-dat=NONE --write-png=ul\.png\$ \


test4d: powerflux.hydra
	valgrind --tool=memcheck --error-limit=no  --leak-check=full ./powerflux.hydra --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--three-bins=0 --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0 --dump-candidates=0 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.7 \
		--write-dat=NONE \
		--write-png=ul\\.png\\$$ 

#./powerflux.hydra --input-format=Power --dataset=test.dst -f $(BAND) -n 501 \
# 		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
# 		--detector=LHO --three-bins=0 --do-cutoff=0 \
# 		--subtract-background=0 --ks-test=1 --filter-lines=0 \
# 		--spindown-start=1e-9 \
# 	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
#         	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat \
# 		--fake-linear --fake-ra=3.141593 --fake-dec=0.006233 --fake-freq=199.85 --fake-strain=1e-20 \
# 		--fake-spindown=1e-9 \

test4e: powerflux
	valgrind --tool=memcheck --error-limit=no  --leak-check=full  ./powerflux --input-format=Power --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --three-bins=0 --do-cutoff=0 \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0e-9 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat \

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

dat2skymap: grid.o rastermagic.o dat2skymap.o png.o
	$(CC) grid.o rastermagic.o dat2skymap.o png.o -o dat2skymap -static $(LDFLAGS)

powerflux: $(OBJS)
	$(CC) $(OBJS) -o powerflux -static $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l powerflux

powerflux.dyn: $(OBJS)
	$(CC) $(OBJS) -o powerflux.dyn $(LDFLAGS) -lgsl -lgslcblas
	ls -l powerflux.dyn

powerflux.hydra: $(OBJS)
	(export PATH=/opt/condor/bin:$$PATH ; \
	 export CONDOR_CONFIG=/opt/condor-6.7.3/etc/condor_config ; \
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
	globus-scp -C powerflux.hydra hydra.phys.uwm.edu:PowerFlux/powerflux.$(VERSION)
	
upload.g: powerflux.hydra
	scp -C powerflux.hydra gallatin.physics.lsa.umich.edu:PowerFlux/powerflux.$(VERSION)
	
upload.h: powerflux.hydra
	globus-scp -C powerflux.hydra ldas-grid.ligo-wa.caltech.edu:PowerFlux/powerflux.$(VERSION)

upload.l: powerflux.hydra
	globus-scp -C powerflux.hydra ldas-grid.ligo-la.caltech.edu:PowerFlux/powerflux.$(VERSION)

upload.c: powerflux.hydra
	globus-scp -C powerflux.hydra ldas-grid.ligo.caltech.edu:PowerFlux/powerflux.$(VERSION)
