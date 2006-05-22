# common rules for application makefiles
# the application makefile needs to define:
# EXEC    = name of the resulting executable
# FEATLIB = list of feat libs to use (feat2d sysutils are required)
# SRC     = list of source files

# optionaly it can define
# INCDIR ... include directive for the compiler ( -Idir )
# DEFS   ... additional defines for compilations (-Dvar )

OBJDIR=obj/$(ID)
MODDIR=obj/$(ID)

#vpath %.inc include
#vpath %.c src
#vpath %.f src

#OBJ =$($(filter %.f,$(SRC)):%.f=$(OBJDIR)/%.o)
#OBJ+=$($(filter %.f90,$(SRC)):%.f90=$(OBJDIR)/%.o)
#OBJ+=$($(filter %.c,$(SRC)):%.c=$(OBJDIR)/%.o)

OBJ=$(filter %.o,$(SRC:%.f=$(OBJDIR)/%.o)) 
OBJ+=$(filter %.o,$(SRC:%.f90=$(OBJDIR)/%.o))
OBJ+=$(filter %.o,$(SRC:%.c=$(OBJDIR)/%.o))

CCOMP=$(CC) $(CCFLAGS) $(OPTFLAGS) $(INCDIR) $(DEFS)
FCOMP=$(FC) $(FCFLAGS) $(OPTFLAGS) $(INCDIR) $(DEFS)
F90COMP=$(FC) $(FCFLAGS) $(OPTFLAGS) $(INCDIR) $(DEFS)


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
	$(F90COMP) $(OBJ) $(LDFLAGS) $(FEAT) $(BLASLIB) -o $@

$(OBJDIR)/%.o : %.f
	$(FCOMP) -c -o $@ $<

$(OBJDIR)/%.o $(MODDIR)/%.mod : %.f90
	$(F90COMP) -c -o $(OBJDIR)/$*.o $<

$(OBJDIR)/%.o : %.c
	$(CCOMP) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

# this checks if the needed libs are ready 
# but it checks only the existence of the lib, if it exists nothing is
# done even if some source files in the lib/src were editited...
$(FEAT):
	@echo "Library" $(@:$(LIBDIR)/lib%.a=%) "is needed by" $(EXEC)
	( cd $(FEATFLOW)/libraries/$(@:$(LIBDIR)/lib%.a=%) && $(MAKE) all )

Deps.mk: $(SRC) $(INC)
	@echo "Recomputing the include dependency file Deps.mk"
	@($(FEATFLOW)/bin/f77mkdep.sh $(filter %.f,$^) >$@)
	@($(FEATFLOW)/bin/f90mkdep.sh $(filter %.f90,$^) >>$@)

clean:
	-rm -f Deps.mk
	-rm -f $(OBJDIR)/*.o 
	-rm -f $(MODDIR)/*.mod 

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

