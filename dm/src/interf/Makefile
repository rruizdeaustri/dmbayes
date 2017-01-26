# Makefile for ge directory
# Author: R. Ruiz


.KEEP_STATE:

.SUFFIXES:      .o .f90


DINC =../../inc
MCMCINC=../../../source

.f90.o:
	$(FC) $(FCFLAGS) -c -I$(DINC)  -I$(MCMCINC) $<

OBJ =  dminit.o dminterface.o intlib.o phig_back.o qshep2d.o sigip.o 

all: $(OBJ)
