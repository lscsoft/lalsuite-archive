PACKAGE=\"powerflux\"
VERSION=1.4.34-64
LALROOT=./../../../../
LALROOT=/opt/lscsoft/lal/
LALINCLUDE=-I$(LALROOT)/include/
LALLDFLAGS=-L$(LALROOT)/lib -L$(LALROOT)/support-lib -L/opt/zlib/lib -L/opt/png/lib -llalsupport -llal -lm -lz
LDFLAGS=$(LALLDFLAGS) -L. -lpng
#OPTFLAGS=-O2 -funroll-loops 
OPTFLAGS=-O6 -mmmx -msse -mcpu=i686
OPTFLAGS=-O6 -march=i686 -funroll-loops
OPTFLAGS=-O6 -msse2 -march=i686 -funroll-loops
OPTFLAGS=-O6 -march=opteron -funroll-loops -fno-stack-protector
OPTFLAGS=-O -march=opteron -fno-stack-protector
OPTFLAGS=-O6 -msse2 -funroll-loops -fno-stack-protector
#OPTFLAGS=-O6 -mtune=opteron -mfpmath=sse -funroll-loops -fno-stack-protector -U_FORTIFY_SOURCE
OPTFLAGS=-O6 -mtune=nocona -funroll-loops -fno-stack-protector -msse3 -U_FORTIFY_SOURCE
CFLAGS=-Wall -g $(OPTFLAGS) -DPACKAGE=$(PACKAGE) -DVERSION=\"$(VERSION)\" \
	$(LALINCLUDE) -fPIC
CC=gcc $(CFLAGS)

COMMON_OBJS=cmdline.o rastermagic.o hookup.o lines.o grid.o \
	png.o polarization.o intervals.o statistics.o dataset.o 

PYTHON_OBJS=$(COMMON_OBJS) candidates.o fine_grid.o

OBJS=$(COMMON_OBJS) fine_grid.o candidates.o powerflux.o

OBJS2=$(COMMON_OBJS) power_cache.o power_sums.o outer_loop.o power_sum_stats.o skymarks.o powerflux2.o summing_context.o

FILES=no_marks.txt sky_marks.txt all_sky_marks_*.txt

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
#DATASET=random.dst
RA=0.4687356
DEC=-1.257051

RA=3.14
DEC=0

#
# Do not build anything by default
#
none:

all: powerflux-condor


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
		--dump-candidates=0 \
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

test4b: powerflux-condor
	valgrind --tool=memcheck --error-limit=no ./powerflux-condor --input-format=Power -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
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


test4d: powerflux-condor
	valgrind --tool=memcheck --error-limit=no  --leak-check=full ./powerflux-condor --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--three-bins=0 --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0 --dump-candidates=0 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.7 \
		--write-dat=NONE \
		--write-png=(Median.*|ul)\\.png\\$$ 

#./powerflux-condor --input-format=Power --dataset=test.dst -f $(BAND) -n 501 \
# 		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
# 		--detector=LHO --three-bins=0 --do-cutoff=0 \
# 		--subtract-background=0 --ks-test=1 --filter-lines=0 \
# 		--spindown-start=1e-9 \
# 	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
#         	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat \
# 		--fake-linear --fake-ra=3.141593 --fake-dec=0.006233 --fake-freq=199.85 --fake-strain=1e-20 \
# 		--fake-spindown=1e-9 \

test4e: powerflux
	valgrind --tool=memcheck --error-limit=no  --leak-check=full  ./powerflux --dataset=$(DATASET) -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LHO --averaging-mode=matched --do-cutoff=0 \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--sky-marks-file=sky_marks.txt \
		--spindown-start=0e-9 \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat

test4f: powerflux.dyn
	( export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LALROOT)/lib:$(LALROOT)/support-lib ; \
	valgrind --tool=memcheck --error-limit=no  --leak-check=full ./powerflux.dyn --dataset=random.dst -f 360000 -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--averaging-mode=one --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=0 --ks-test=1 --filter-lines=0 \
		--spindown-start=0.0 --spindown-step=2e-10 --spindown-count=2 \
		--dump-candidates=0 \
		--dump-points=0 \
		--sky-marks-file=sky_marks.txt \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
 		--write-dat=NONE \
 		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.1 \
 		--fake-ref-time=793154935 \
 		--fake-ra=$(RA) --fake-dec=$(DEC) --fake-freq=200.10027777 --fake-strain=10e-24 \
 		--fake-spindown=0e-10 --fake-iota=1.57 --fake-psi=0.0 \
		--max-candidates=10000 --skymap-resolution-ratio=1 \
	)

test4g: powerflux-mt
	( export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LALROOT)/lib:$(LALROOT)/support-lib ; \
	valgrind --tool=memcheck --error-limit=no  --leak-check=full ./powerflux-mt --dataset=random.dst -f 360000 -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--averaging-mode=one --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=0 --ks-test=1 --filter-lines=0 \
		--spindown-start=0.0 --spindown-step=2e-10 --spindown-count=2 \
		--dump-candidates=0 \
		--dump-points=0 \
		--sky-marks-file=all_sky_marks.txt \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
 		--write-dat=NONE \
 		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.1 \
		--max-candidates=10000 --skymap-resolution-ratio=1 \
	)

test5: powerflux
	time ./powerflux -i /home/volodya/LIGO/S2/S2.L1.op/power.L1. -f $(BAND) -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--detector=LLO --three-bins=0 --do-cutoff=0 --filter-lines=1 \
		--fake-strain=0.3e-22 --config=test_conf --resolution-ratio=1.0 \
		--fake-orientation=0.5 --fake-ra=5.785685 --fake-dec=-1.563849

test5c: powerflux-mt
	 ./powerflux-mt --dataset=random.dst -f 180002 -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--averaging-mode=single --do-cutoff=1 --lock-file=test.lock \
		--subtract-background=1 --ks-test=1 --filter-lines=0 \
		--spindown-start=0.0 --spindown-step=2e-8 --spindown-count=1 \
		--dump-candidates=0 \
		--dump-points=0 --no-secondary-skymaps=1 \
		--sky-marks-file=all_sky_marks.txt \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--max-candidates=10000 --skymap-resolution-ratio=1 \
  		--focus-dec=$(DEC) --focus-ra=$(RA) --focus-radius=0.5e2 

#  		--fake-ref-time=793154935 \
#  		--fake-ra=$(RA) --fake-dec=$(DEC) --fake-freq=200.10027777 --fake-strain=10e-24 \
#  		--fake-spindown=2e-10 --fake-iota=1.57 --fake-psi=0.0 \


# 		--write-dat=NONE \
#		--write-png=NONE \
# 		--write-png=\(.\*bands\|.\*strain.\*\|high.\*\|S.\*\|Median.\*\|weight.\*\|max_dx.\*\)\\.png$$ \
#		--dump-data=input.dat \

test_bench_one: powerflux-mt
	./powerflux-mt --dataset=random3.dst -f 180005 -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--averaging-mode=one --do-cutoff=1 \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0.0 --spindown-step=2e-8 --spindown-count=1 \
		--dump-candidates=0 \
		--dump-points=0 \
		--sky-marks-file=all_sky_marks.txt \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--max-candidates=10000 --skymap-resolution-ratio=1 \
		--num-threads=2

test_bench_matched: powerflux-mt
	./powerflux-mt --dataset=random3.dst -f 180007 -n 501 \
		--ephemeris-path=/home/volodya/LIGO/LAL/lalapps/src/detresponse \
		--averaging-mode=matched --do-cutoff=1 \
		--subtract-background=1 --ks-test=1 --filter-lines=1 \
		--spindown-start=0.0 --spindown-step=2e-8 --spindown-count=1 \
		--dump-candidates=0 \
		--dump-points=0 \
		--sky-marks-file=all_sky_marks.txt \
	        --earth-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/earth05-09.dat \
        	--sun-ephemeris=/home/volodya/LIGO/LAL/lal/packages/pulsar/test/sun05-09.dat  \
		--max-candidates=10000 --skymap-resolution-ratio=1 \
		--num-threads=2


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

python.o: python.c *.h
	$(CC) -c `python-config --includes` python.c

powerflux.so: $(PYTHON_OBJS) jobs-mt.o util.o python.o
	$(CC) -shared $(PYTHON_OBJS) jobs-mt.o util.o python.o -o powerflux.so $(LALLDFLAGS) /usr/lib/libpng.so -lgsl -lgslcblas -lm
	ls -l powerflux.so

powerflux: $(OBJS) jobs.o util.o
	$(CC) $(OBJS) jobs.o util.o -o powerflux -static $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l powerflux

powerflux-mt: $(OBJS) jobs-mt.o util.o
	$(CC) $(OBJS) jobs-mt.o util.o -o powerflux-mt -static $(LDFLAGS) -lgsl -lgslcblas -lm -lpthread -lz
	ls -l powerflux-mt

powerflux.dyn: $(OBJS) jobs.o util.o
	$(CC) $(OBJS) jobs.o util.o -o powerflux.dyn $(LDFLAGS) -lgsl -lgslcblas -lz
	ls -l powerflux.dyn

powerflux-condor: $(OBJS) jobs.o util-condor.o
	 condor_compile $(CC) $(OBJS) jobs.o util-condor.o -o powerflux-condor -static $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l powerflux-condor

powerflux-mt2: $(OBJS2) jobs-mt.o util.o
	$(CC) $(OBJS2) jobs-mt.o util.o -o powerflux-mt2 -static $(LDFLAGS) -lgsl -lgslcblas -lm -lpthread -lz
	ls -l powerflux-mt2

powerflux-condor2: $(OBJS2) jobs.o util-condor.o
	condor_compile $(CC) $(OBJS2) jobs.o util-condor.o -o powerflux-condor2 -static $(LDFLAGS) -lgsl -lgslcblas -lm
	ls -l powerflux-condor2

powerflux-valgrind2: $(OBJS2) jobs-mt.o util.o
	$(CC) $(OBJS2) jobs-mt.o util.o -o powerflux-valgrind2 $(LDFLAGS) -lgsl -lgslcblas -lm -lpthread -lz
	ls -l powerflux-valgrind2


jobs-mt.o: jobs.c jobs.h
	$(CC) $(CFLAGS) -DLINUX_THREADS -c jobs.c -o jobs-mt.o

util-condor.o: util.c util.h
	$(CC) $(CFLAGS) -DCONDOR -c util.c -o util-condor.o

test_rm: test_rm.c rastermagic.o
	$(CC) test_rm.c rastermagic.o -lm -o test_rm

cmdline.c: powerflux.ggo
	gengetopt  --string-parser --conf-parser < powerflux.ggo

cmdline2.c: powerflux2.ggo
	gengetopt  --string-parser --conf-parser --file-name=cmdline2 < powerflux2.ggo

clean:
	rm -f core *.bck *.o

upload: powerflux-condor powerflux-mt
	globus-scp -C powerflux-condor hydra.phys.uwm.edu:PowerFlux/powerflux-condor.$(VERSION)
	globus-scp -C powerflux-mt hydra.phys.uwm.edu:PowerFlux/powerflux-mt.$(VERSION)
	globus-scp -C $(FILES) hydra.phys.uwm.edu:PowerFlux/
	
upload.a: powerflux-condor powerflux-mt powerflux-condor2
	globus-scp -C powerflux-condor h2.atlas.aei.uni-hannover.de:PowerFlux/powerflux-condor.$(VERSION)
	globus-scp -C powerflux-mt h2.atlas.aei.uni-hannover.de:PowerFlux/powerflux-mt.$(VERSION)
	globus-scp -C powerflux-condor2 h2.atlas.aei.uni-hannover.de:PowerFlux/powerflux-condor2.$(VERSION)
	globus-scp -C $(FILES) h2.atlas.aei.uni-hannover.de:PowerFlux/
	
upload.g: powerflux-condor powerflux-mt
	scp -C powerflux-condor gallatin.physics.lsa.umich.edu:PowerFlux/powerflux-condor.$(VERSION)
	scp -C powerflux-mt gallatin.physics.lsa.umich.edu:PowerFlux/powerflux-mt.$(VERSION)
	scp -C $(FILES) gallatin.physics.lsa.umich.edu:PowerFlux/
	
upload.h: powerflux-condor powerflux-mt
	globus-scp -C powerflux-condor ldas-grid.ligo-wa.caltech.edu:PowerFlux/powerflux-condor.$(VERSION)
	globus-scp -C $(FILES) ldas-grid.ligo-wa.caltech.edu:PowerFlux/

upload.l: powerflux-condor powerflux-mt
	globus-scp -C powerflux-condor ldas-grid.ligo-la.caltech.edu:PowerFlux/powerflux-condor.$(VERSION)
	globus-scp -C $(FILES) ldas-grid.ligo-la.caltech.edu:PowerFlux/

upload.c: powerflux-condor powerflux-mt powerflux-condor2
	globus-scp -C powerflux-condor ldas-grid.ligo.caltech.edu:PowerFlux/powerflux-condor.$(VERSION)
	globus-scp -C powerflux-mt ldas-grid.ligo.caltech.edu:PowerFlux/powerflux-mt.$(VERSION)
	globus-scp -C powerflux-condor2 ldas-grid.ligo.caltech.edu:PowerFlux/powerflux-condor2.$(VERSION)
	globus-scp -C $(FILES) ldas-grid.ligo.caltech.edu:PowerFlux/

upload.m: powerflux-condor powerflux-mt
	globus-scp -C powerflux-condor vdergach@morgane.aei.mpg.de:PowerFlux/powerflux-condor.$(VERSION)
	globus-scp -C powerflux-mt vdergach@morgane.aei.mpg.de:PowerFlux/powerflux-mt.$(VERSION)
	globus-scp -C $(FILES) vdergach@morgane.aei.mpg.de:PowerFlux/
	
