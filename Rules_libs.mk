# common rules for libraries makefiles
# the library makefile needs to define: 
# LIBNAME = name of the resulting library
# SRCLIST = list of source files
#
# optionally it can define
# INCDIR  = include directive for the compiler ( -Idir )
# DEFS    = additional defines for compilations (-Dvar )


DIR=$(FEATFLOW)/$(LIBNAME)
LIB=$(LIBNAME:%=$(LIBDIR)/lib%.a)

SRCDIR=src
OBJDIR=obj/$(ID)

#vpath %.inc include
#vpath %.c src
#vpath %.f src

SRC=$(SRCLIST:%=$(SRCDIR)/%)
OBJ=$(filter %.o,$(SRCLIST:%.f=$(OBJDIR)/%.o)) 
OBJ+=$(filter %.o,$(SRCLIST:%.c=$(OBJDIR)/%.o))

CCOMP=$(CC) $(CCFLAGS) $(OPTFLAGS) $(INCDIR) $(DEFS)
FCOMP=$(FC) $(FCFLAGS) $(OPTFLAGS) $(INCDIR) $(DEFS)

all: greet lib
	@echo "Done," $(LIBNAME) "is ready."
	@echo

lib: $(OBJDIR) $(LIB)

greet:
	@echo "Compiling library" $(LIBNAME) "in" 
	@pwd

$(LIB): $(OBJ)
	mkdir -p $(LIBDIR)
	$(AR) -sr $@ $? 
	ranlib $@

$(OBJDIR)/%.o : $(SRCDIR)/%.f
	$(FCOMP) -c -o $@ $<

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CCOMP) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	-rm -f $(OBJDIR)/*.o 
	-rm -f Deps.mk

purge: clean
	-rm -f $(LIB)

purge_all: clean
	-rm -f -r obj/*
	-rm -f $(FEATFLOW)/object/libraries/lib-*/$(LIBNAME:%=lib%.a)

id: .id
help: .help
debug: all

