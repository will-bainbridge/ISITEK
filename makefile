COMPILER	= gcc
FLAG		= -g -Wall

HOMEPATH	= $(shell pwd)
BUILDPATH	= $(HOMEPATH)/src
THIRDPATH	= $(HOMEPATH)/thirdparty

#------------------------------------------------------------------------------#

INCLUDE		+= -I$(THIRDPATH)/UMFPACK/Include -I$(THIRDPATH)/AMD/Include -I$(THIRDPATH)/UFconfig
LIBRARY		+= -L$(THIRDPATH)/UMFPACK/Lib -L$(THIRDPATH)/AMD/Lib -lumfpack -lamd

LIBRARY		+= -llapack -lblas
#LIBRARY		+= -L$(THIRDPATH)/GotoBLAS2 -Wl,-R$(THIRDPATH)/GotoBLAS2 -lgoto2

LIBRARY		+= -lm -lrt -lgfortran

#------------------------------------------------------------------------------#

MAINFILES	= isitek.c
COMMONFILES	= expression.c fetch.c geometry.c io.c memory.c numerics.c sparse.c system.c

MAINSOURCE	= $(addprefix $(BUILDPATH)/,$(MAINFILES))
COMMONSOURCE 	= $(addprefix $(BUILDPATH)/,$(COMMONFILES))
ALLSOURCE	= $(MAINSOURCE) $(COMMONSOURCE)

MAINOBJECT	= $(MAINSOURCE:.c=.o)
COMMONOBJECT	= $(COMMONSOURCE:.c=.o)
ALLOBJECT	= $(MAINOBJECT) $(COMMONOBJECT)

EXECUTABLES	= $(MAINFILES:.c=)

################################################################################

all: $(EXECUTABLES)

.SECONDEXPANSION:
$(EXECUTABLES): $(BUILDPATH)/$$@.o $(COMMONOBJECT)
	$(COMPILER) $(FLAG) -o $@ $(BUILDPATH)/$@.o $(COMMONOBJECT) $(LIBRARY)

$(ALLOBJECT): makefile
	$(COMPILER) $(FLAG) -c $*.c -o $*.o $(INCLUDE)

-include $(BUILDPATH)/depend
depend: $(ALLSOURCE)
	$(COMPILER) -MM $(INCLUDE) $^ | sed 's|^\(.*\.o\)|$(BUILDPATH)/\1|g' > $(BUILDPATH)/$@

clean:
	rm -f $(BUILDPATH)/*o
