
FRAMEROOT=/home/volodya/LIGO/framelib/
LALROOT=./../../../../../
LALINCLUDE=-I$(LALROOT)/include/
LALLDFLAGS=-L$(LALROOT)/lib -L$(LALROOT)/support-lib -llalsupport -llal -lm


OPTFLAGS=-O2 -mcpu=i686 -funroll-loops

CFLAGS=-Wall -g $(OPTFLAGS) $(LALINCLUDE) -I$(FRAMEROOT)/include/

LDFLAGS=$(LALLDFLAGS) -L$(FRAMEROOT)/bin/ 


#CFLAGS=-Wall -g -O2 -I../framelib/include -I../LAL/include
#LDFLAGS=-L../framelib/bin/ -L../LAL/lib/ -L/ldcg/lib -L../LAL/support-lib/

CC=gcc -pipe $(CFLAGS)

UPLOAD_FILES1=README.txt Makefile make_sft.c weighted_sum.c make_sft_conf.l \
	make_sft_global.h make_sft_calibration.c input.example \
	make_cross_section.c cross_section.tcl op_method.c op_method.h \
	make_sft_utils.c make_sft_td_data.c make_sft_op.c make_sft_op_conf.l
	
UPLOAD_FILES=make_sft_op make_sft	
	
TAR_FILES=SFT-3/README.txt SFT-3/Makefile SFT-3/make_sft.c SFT-3/weighted_sum.c \
	 SFT-3/make_sft_conf.l SFT-3/make_sft_global.h SFT-3/input.example \
	 SFT-3/make_sft_calibration.c SFT-3/make_cross_section.c SFT-3/cross_section.tcl \
	 SFT-3/op_method.c SFT-3/op_method.h SFT-3/make_sft_utils.c \
	 SFT-3/make_sft_td_data.c SFT-3/make_sft_op.c SFT-3/make_sft_op_conf.l

all: 
# Do nothing so lalapps build does not break

programs: make_sft weighted_sum make_cross_section make_sft_op


example: make_sft
	time ./make_sft < input.example

test1: make_sft
	time ./make_sft < input.test1
	
test2: make_sft
	time ./make_sft < input.test2
	
testop1: make_sft_op
	time ./make_sft_op < input.testop1
	
testw: weighted_sum
	time ./weighted_sum 0.333333 out-1/out.1.gz 0.3333333 out-1/out.2.gz 0.333333 out-2/out.LLO.1 | less 

test_ag: adjust_geo
	cp sft.H1.1 sft.H1.1t
	./adjust_geo sft.H1.1t

adjust_geo: adjust_geo.c
	$(CC) $(CFLAGS) $(LDFLAGS) adjust_geo.c -o adjust_geo

make_sft_conf.c: make_sft_conf.l
	flex -omake_sft_conf.c make_sft_conf.l
	
make_sft_op_conf.c: make_sft_op_conf.l
	flex -omake_sft_op_conf.c make_sft_op_conf.l

	
make_sft: make_sft.c make_sft_conf.c make_sft_global.h make_sft_calibration.c op_method.h op_method.c\
	make_sft_utils.c make_sft_td_data.c
	$(CC) $(CFLAGS) $(LDFLAGS) make_sft.c make_sft_conf.c make_sft_calibration.c \
	op_method.c make_sft_utils.c make_sft_td_data.c \
	-o make_sft -static -lFrame -llal -lm -lfftw3f -lfftw3
	ls -l make_sft

make_sft_op: make_sft_op.c make_sft_op_conf.c make_sft_global.h make_sft_calibration.c op_method.h op_method.c\
	make_sft_utils.c make_sft_td_data.c
	$(CC) -static $(CFLAGS) $(LDFLAGS) make_sft_op.c make_sft_op_conf.c make_sft_calibration.c \
	op_method.c make_sft_utils.c make_sft_td_data.c \
	-o make_sft_op  -lFrame  -llal -llalsupport -lm  -lfftw3f -lfftw3
	ls -l make_sft_op

weighted_sum: weighted_sum.c
	$(CC) $(CFLAGS) $(LDFLAGS) weighted_sum.c -lm -o weighted_sum

make_cross_section: make_cross_section.c
	$(CC) $(CFLAGS) $(LDFLAGS) make_cross_section.c -lm -o make_cross_section

condor: make_sft_condor

make_sft_condor: make_sft.c make_sft_conf.c make_sft_global.h make_sft_calibration.c
	condor_compile gcc -Wall -pipe -O -I/ldcg/include -L/ldcg/lib -L../LAL/lib -I../LAL/include make_sft.c make_sft_conf.c make_sft_calibration.c -o make_sft_condor -static -lFrame -llal -lsrfftw -lsfftw  -lm

hydra: make_sft_hydra make_sft_op_hydra

make_sft_hydra: make_sft.c make_sft_conf.c make_sft_global.h make_sft_calibration.c \
	op_method.h op_method.c make_sft_utils.c
	gcc $(CFLAGS) -L../framelib/bin -L../LAL/lib  make_sft.c make_sft_conf.c make_sft_calibration.c \
	op_method.c make_sft_utils.c make_sft_td_data.c \
	-o make_sft_hydra -static -lFrame -llal -lsrfftw -lsfftw -lm


make_sft_op_hydra: make_sft_op.c make_sft_op_conf.c make_sft_global.h make_sft_calibration.c op_method.h op_method.c\
	make_sft_utils.c make_sft_td_data.c
	$(CC) $(CFLAGS) -L../framelib/bin -L../LAL/lib make_sft_op.c make_sft_op_conf.c make_sft_calibration.c \
	op_method.c make_sft_utils.c make_sft_td_data.c \
	-o make_sft_op_hydra -static -lFrame -llal  -lsrfftw -lsfftw  -lm 

upload: $(UPLOAD_FILES)
	scp -C -S globus-ssh $(UPLOAD_FILES) hydra.phys.uwm.edu:/home/volodya/SFT-3/

upload.h: $(UPLOAD_FILES)
	scp -C -S globus-ssh $(UPLOAD_FILES) ldas-grid.ligo-wa.caltech.edu:SFT-3/

upload-tenaya: $(UPLOAD_FILES)
	scp -C $(UPLOAD_FILES) tenaya.physics.lsa.umich.edu:/home/volodya/SFT-3/

tarball:
	(cd .. ; tar cvf - $(TAR_FILES) ) | gzip - > SFT-3.tgz
