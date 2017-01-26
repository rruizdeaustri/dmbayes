# Makefile template for Diver 1.0.0 beta 
# Author: Pat Scott, patscott@physics.mcgill.ca

# Until we sort out some configure magic, change 
# the options below to suit your system

# Define fortran compiler and options
FF=mpif90
FOPT=-DMPI -O3 -fPIC -ffree-line-length-none -ffixed-line-length-none -cpp

export FF FOPT

DIVE_ROOT = $(PWD)

AR = ar r  
LINKLIB = ld -shared  
SOURCEFILES = detypes deutils mutation crossover selection init converge posterior evidence io de cwrapper
SOURCE = ${DIVE_ROOT}/source
INC = ${DIVE_ROOT}/include
BUILD = ${DIVE_ROOT}/build

OBJ = $(SOURCEFILES:%=$(BUILD)/%.o)

all: libdiver.a

libdiver.a: $(OBJ) makefile
	$(AR) $@ $(OBJ) 

libdiver.so: $(OBJ) makefile
	$(LINKLIB) -o $@ $(OBJ) 

$(BUILD)/%.o: $(SOURCE)/%.f90
	cd $(BUILD); \
	$(FF) -c $(FOPT) -I$(INC) $<

$(BUILD)/%.o: $(SOURCE)/%.F90
	cd $(BUILD); \
	$(FF) -c $(FOPT) -I$(INC) $<

clean:
	rm -f *.a *.so; \
	cd $(BUILD); rm -f *.o *.mod 
