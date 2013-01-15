################################################################################

COMPILER	= gcc
FLAG		= -g -Wall

HOMEPATH	= .
BUILDPATH	= $(HOMEPATH)/src
THIRDPATH	= $(HOMEPATH)/thirdparty
THIRDFULLPATH	= $(shell cd $(THIRDPATH);pwd)

#------------------------------------------------------------------------------#

MAINFILES	= isitek.c
COMMONFILES	= boundary.c condition.c element.c element_interpolation.c expression.c face.c face_interpolation.c memory.c node.c numerics.c quadrature.c solver.c sparse.c

MAINSOURCE	= $(addprefix $(BUILDPATH)/,$(MAINFILES))
COMMONSOURCE 	= $(addprefix $(BUILDPATH)/,$(COMMONFILES))
ALLSOURCE	= $(MAINSOURCE) $(COMMONSOURCE)

MAINOBJECT	= $(MAINSOURCE:.c=.o)
COMMONOBJECT	= $(COMMONSOURCE:.c=.o)
ALLOBJECT	= $(MAINOBJECT) $(COMMONOBJECT)

EXECUTABLES	= $(MAINFILES:.c=)

#------------------------------------------------------------------------------#

LIBRARY += -lm -lrt

#LIBRARY += -lblas -llapack

$(BUILDPATH)/sparse.o depend: FLAG += -DSOLVE_UMFPACK
$(BUILDPATH)/sparse.o depend: INCLUDE += -I$(THIRDPATH)/UMFPACK/Include -I$(THIRDPATH)/AMD/Include -I$(THIRDPATH)/UFconfig
LIBRARY += -L$(THIRDPATH)/UMFPACK/Lib -L$(THIRDPATH)/AMD/Lib -lumfpack -lamd
LIBRARY += -L$(THIRDPATH)/GotoBLAS2 -Wl,-R$(THIRDFULLPATH)/GotoBLAS2 -lgoto2 -lgfortran

#$(BUILDPATH)/sparse.o depend: FLAG += -DSOLVE_PARDISO
#$(BUILDPATH)/sparse.o depend: INCLUDE += -I/opt/intel/mkl/include
##LIBRARY += -L/opt/intel/mkl/lib/intel64 -Wl,-R/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_solver_lp64
##LIBRARY += -L/opt/intel/composerxe-2011.1.107/compiler/lib/intel64 -liomp5
#LIBRARY += -L/opt/intel/mkl/lib/intel64 -Wl,-R/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_solver_lp64
#LIBRARY += -lpthread

################################################################################

all: $(EXECUTABLES)

.SECONDEXPANSION:
$(EXECUTABLES): $(BUILDPATH)/$$@.o $(COMMONOBJECT)
	$(COMPILER) $(FLAG) -o $@ $(BUILDPATH)/$@.o $(COMMONOBJECT) $(LIBRARY)

$(ALLOBJECT): makefile
	$(COMPILER) $(FLAG) -c $*.c -o $*.o $(INCLUDE)

-include $(BUILDPATH)/depend
depend: $(ALLSOURCE)
	$(COMPILER) $(FLAG) -MM $(INCLUDE) $^ | sed 's|^\(.*\.o\)|$(BUILDPATH)/\1|g' > $(BUILDPATH)/$@

clean:
	rm -f $(BUILDPATH)/*o

################################################################################
