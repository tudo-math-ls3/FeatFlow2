# common rules for application makefiles
# the application makefile needs to define:
# EXEC    = name of the resulting executable
# FEATLIB = list of feat libs to use (feat2d sysutils are required)
# SRC     = list of source files


OBJDIR=obj/$(ID)
OBJ=$(SRC:%.f=$(OBJDIR)/%.o)

#vpath %.inc include
#vpath %.c src
#vpath %.f src

# If the BLASLIB is not defined add the included blas to the libs.
# If the LAPACKLIB is not defined add the included lapack to the libs.
FEATLIB+= $(if $(LAPACKLIB),,lapack) $(if $(BLASLIB),,blas)

FEAT=$(FEATLIB:%=$(LIBDIR)/lib%.a)

all: greet $(EXEC)
	@echo "Done," $(EXEC) "is ready."
	@echo

greet:
	@echo "Compiling module" $(EXEC) "in" 
	@pwd

# If a BLAS or LAPACK library is given, BLASLAPACKLIB contains entries
# for both, BLAS and LAPACK.

$(EXEC): $(OBJDIR) $(FEAT) Deps.mk $(OBJ) 
	$(LD) $(FCFLAGS) $(OPTFLAGS) $(OBJ) $(LDFLAGS) $(FEAT) $(BLASLAPACKLIB) $(LDLIBS) -o $@

# The variable BENCHINC is probably set by the benchmark if the user
# executes it. This variable adds additionally INCLUDE directories
# to the call of the compiler, otherwise some INC files would not
# be found!

$(OBJDIR)/%.o : %.f
	$(FC) $(FCFLAGS) $(OPTFLAGS) $(INCDIR) $(BENCHINC) $(DEFS) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

# this checks if the needed libs are ready 
# but it checks only the existence of the lib, if it exists nothing is
# done even if some source files in the lib/src were editited...
$(FEAT):
	@echo "Library" $(@:$(LIBDIR)/lib%.a=%) "is needed by" $(EXEC)
	( cd $(FEATFLOW)/libraries/$(@:$(LIBDIR)/lib%.a=%) && $(MAKE) all )

Deps.mk: $(SRC)
	($(FEATFLOW)/bin/f77mkdep.sh $^ >$@)

clean:
	-rm -f Deps.mk
	-rm -f $(OBJDIR)/*.o 

clean_exec: 
	-rm -f $(EXEC)

purge: clean clean_exec

purge_all: purge
	@(for i in *-*-*-* ; do if [ -x $$i ] ; then echo rm $$i ; rm -f $$i ; fi ; done )
	-rm -f -r obj/*
	-rm -f -r \#gmv/* \#avs/* \#points/* \#ns/* \#film/* \#data/\#DX*
	-rm -f -r \#data/*.ERR \#data/*.PRT \#data/*.SYS \#data/*.TRC 

.NOTPARALLEL: clean purge run

debug: all

run: all
	@(time ./$(EXEC))

id: .id
help: .help
-include Deps.mk

