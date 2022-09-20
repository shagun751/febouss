SWITCH_ICC_RELEASE=-ipo -O3 -vec_report1
SWITCH_ICC_CUSTOM=-O3
SWITCH_ICC_DEBUG= -g -debug all -check all -implicitnone -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -traceback
SWITCH_GFORTRAN_RELEASE=-Ofast
SWITCH_GFORTRAN_O2=-O2
SWITCH_GFORTRAN_DEBUG=-O0 -g -fbounds-check -fbacktrace

#COMPILER=ifort
COMPILER=gfortran


ifeq ($(COMPILER),ifort)
	MODPATH=-module subroutines
	C_COMP=icc
	PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64
	PARA_PATH=$(HOME)/Programs/paralution_icc_min
else ifeq ($(COMPILER),gfortran)
	MODPATH=-J subroutines
	C_COMP=gcc
	PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64
	PARA_PATH=$(HOME)/Programs/paralution_gcc_min
endif
OTH_FLAGS = -cpp 

PARA_CONV=./subroutines/paralution_fortran.cpp
PARA_LIB=$(PARA_PATH)/build/lib/libparalution.so 
#PARA_LIB=./libparalution.so 
PARA_INC=-I$(PARA_PATH)/build/inc

MAIN=bsnqQuadMain
MODULE=bsnqModule

OBJ=paralution_fortran.o $(MODULE).o $(MAIN).f90

SUBS=\
subroutines/$(MODULE).f90 \
subroutines/bsnqModuleFncs.f90 \
subroutines/bsnqModuleFncs2.f90 \
subroutines/initialCondition.f90 \
subroutines/outputXML.f90 \
subroutines/mods.f90 \
subroutines/modsInletBC.f90 \
subroutines/modsMFree.f90 \
subroutines/modsVertVel.f90 \
subroutines/shipMod.f90 \
subroutines/bndCondition_v1.0.f90 \
subroutines/bndIntegral.f90 \
subroutines/femAnalyticalTri_v1.1.f90 \
subroutines/findNeiFromLinkList.f90 \
subroutines/geometry_v1.0.f90 \
subroutines/matrixSet1.f90 \
subroutines/matrixSet2.f90 \
subroutines/mergeSort.f90 \
subroutines/nodeConnAll.f90 \
subroutines/rk4Interpolation.f90 \
subroutines/solver_v1.0.f90 \
subroutines/modCellStore.f90
	
default: $(OBJ)
	mkdir -p Output
	$(COMPILER) -o bsnq $(OBJ) $(PARA_LIB) $(PARA_FLAGS) $(MODPATH) $(OTH_FLAGS)
	
# $(MAIN).o: $(MAIN).f90
# 	ifort $(PARA_FLAGS) -c $(MAIN).f90 $(SUBS)

#g++ $(PARA_FLAGS) $(PARA_INC) -c $(PARA_CONV)

$(MODULE).o: $(SUBS)
	$(COMPILER) -c $< $(PARA_FLAGS) $(MODPATH)

paralution_fortran.o: $(PARA_CONV)	
	$(C_COMP) $(PARA_FLAGS) $(PARA_INC) -c $(PARA_CONV)	

clean:
	rm -rf Output 
	mkdir Output

cleanAll:
	rm -rf *.o subroutines/*.mod *.mod Output bsnq 
	mkdir Output

rmout:
	rm -rf Output
	mkdir Output
	
run: clean default	
	./bsnq rect2D
	cat rect2D.rout


