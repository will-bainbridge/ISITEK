COMPILER	= gcc
FLAG		= -O2 -Wall

HOMEPATH	= .
BUILDPATH	= $(HOMEPATH)/src
THIRDPATH	= $(HOMEPATH)/thirdparty
THIRDFULLPATH	= $(shell cd $(THIRDPATH);pwd)

#------------------------------------------------------------------------------#

MAINFILES	= isitek.c
COMMONFILES	= constants.c expression.c fetch.c geometry.c io.c memory.c numerics.c sparse.c system.c

MAINSOURCE	= $(addprefix $(BUILDPATH)/,$(MAINFILES))
COMMONSOURCE 	= $(addprefix $(BUILDPATH)/,$(COMMONFILES))
ALLSOURCE	= $(MAINSOURCE) $(COMMONSOURCE)

MAINOBJECT	= $(MAINSOURCE:.c=.o)
COMMONOBJECT	= $(COMMONSOURCE:.c=.o)
ALLOBJECT	= $(MAINOBJECT) $(COMMONOBJECT)

EXECUTABLES	= $(MAINFILES:.c=)

#------------------------------------------------------------------------------#

LIBRARY += -lm -lrt

$(BUILDPATH)/sparse.o depend: FLAG += -DSOLVE_UMFPACK
$(BUILDPATH)/sparse.o depend: INCLUDE += -I$(THIRDPATH)/SuiteSparse/SuiteSparse_config -I$(THIRDPATH)/SuiteSparse/AMD/Include -I$(THIRDPATH)/SuiteSparse/UMFPACK/Include
LIBRARY += -L$(THIRDPATH)/SuiteSparse/lib -Wl,-R$(THIRDFULLPATH)/SuiteSparse/lib -lmetis -lsuitesparseconfig -lamd -lcamd -lcolamd -lccolamd -lcholmod -lumfpack
LIBRARY += -L$(THIRDPATH)/OpenBLAS -Wl,-R$(THIRDFULLPATH)/OpenBLAS -lopenblas -lgfortran

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
