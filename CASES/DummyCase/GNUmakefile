#AMREX_HOME    = $(HOME)/Softwares/amrex/build2d
#HYPRE_DIR     = /usr/local/hypre
#BLITZ_DIR     = $(HOME)/software/blitz
#MPI_DIR       = /usr/lib/x86_64-linux-gnu/openmpi
#ALGOIM        = $(HOME)/Softwares/algoim/src


AMREX_HOME      = /home/pratik/NumericalLibraries/amrex-latest/amrex-development/tmp_install_dir
#AMREX_HOME      = /home/pratik/NumericalLibraries/amrex/tmp_install_dir
HYPRE_DIR      = /home/pratik/NumericalLibraries/hypre/src/hypre
ALGOIM        = /home/pratik/NumericalLibraries/algoim/src
MPI_DIR       = /usr/bin/
BLITZ_DIR     = /home/pratik/NumericalLibraries/blitz_installation/
PETSC_DIR     = /home/pratik/NumericalLibraries/petsc-3.7.7

SRC = ../../Src
APP = ../../Application/incFSIAMR

##INC_PATH      = -I$(AMREX_HOME)/include -I$(HYPRE_DIR)/include -I$(SRC) -I$(SRC)/lninclude \
                -I$(ALGOIM) -I$(PETSC_DIR)/include -I$(APP) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/gnu-opt/include -I$(MPI_DIR)/include
# -I$(MPI_DIR)/include
INC_PATH      = -I$(AMREX_HOME)/include -I$(HYPRE_DIR)/include -I$(SRC) -I$(SRC)/lninclude \
                -I$(ALGOIM) -I$(APP) -I$(MPI_DIR)/include

LIB_PATH      = -L$(AMREX_HOME)/lib -L$(HYPRE_DIR)/lib 

COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(AMREX_HOME)/lib/pkgconfig/amrex.pc)
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(AMREX_HOME)/lib/pkgconfig/amrex.pc)

CFLAGS := $(INC_PATH) $(COMPILE_CPP_FLAGS)
LFLAGS := $(COMPILE_LIB_FLAGS)

CXX = mpicxx
CXXFLAGS = -O2 -std=c++14 -g -Wall

#csrc = $(wildcard $(foreach dir, $(Pdirs), $(SRC)/$(dir)/*.cpp) )
csrc += $(wildcard $(SRC)/**/*.cpp)
csrc += $(wildcard $(APP)/*.cpp)
csrc += $(wildcard ./*.cpp)
obj = $(csrc:.cpp=.o)
dep = $(obj:.o=.d)

default: run.exe

run.exe: $(obj)
	$(CXX) -o $@ $(CXXFLAGS) $(CFLAGS) $^ $(LIB_PATH) $(LFLAGS) -lHYPRE 
#-lblitz
#-lspdlog

%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $(CFLAGS) $<

-include $(dep)

%.d: %.cpp
	@$(CPP) $(CXXFLAGS) $(INC_PATH) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean realclean

clean:
	$(RM) $(SRC)/*.o

realclean:
	$(RM) $(obj) $(dep)
	$(RM) run.exe

all:default


srun:
	./run.exe inputs
##parrun:
##	/home/pratik/NumericalLibraries/petsc-3.7.7/gnu-opt/bin/mpirun -np 4 ./run.exe inputs
parrun:
	mpirun -np 6 ./run.exe inputs
valsrun:
	valgrind --log-file="val.out" ./run.exe inputs
cleandata:
	rm -rf Output/*
