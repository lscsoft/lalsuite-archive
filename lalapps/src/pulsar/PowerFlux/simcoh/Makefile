#
#    Template make file to link DMT monitors on a GC computer. All libraries needed for
#    the DMT data environment (DatEnv base class) are included in the 
#    loader library list. 
#
SHELL		:= /bin/sh

#--->  Default version.
ifndef DMTVERSION
DMTVERSION	:= pro
endif

INC		:= $(DMTHOME)/$(DMTVERSION)/include
LIB		:= $(DMTHOME)/$(DMTVERSION)/lib

#CC		:= /export/apps4/ldcg/bin/gcc
#CXX		:= /export/apps4/ldcg/bin/g++
#LD		:= /export/apps4/ldcg/bin/g++
CC              := gcc
CXX             := g++
LD              := g++
CFLAGS	:= -Wall -gstabs -I/home/jyzhang/myinclude -I /home/rarmen/opt/lscsoft/non-lsc/include/ -I$(LAL_LOCATION)/include
LDLIBS          := -llalsupport -llal -lm
#
#   monitors using the root libraries should add the following: 
#



#
#   Monitors linking together multiple object files should list the 
#   dependencies as follows:
#
simcoh: 	simcoh.o 

#
#
clean:
	rm *.o simcoh.o

#
#   and specify programs to be linked in 'PROGRAMS'
#
PROGRAMS	:= simcoh
#
ifdef PROGRAMS
$(PROGRAMS):	
	$(LD) $^ -o $@ $(LDLIBS)
endif
