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
OBJ+=$(filter %.o,$(SRC:%.cpp=$(OBJDIR)/%.o))

CCOMP=$(CC) $(CCFLAGS) $(OPTFLAGS) $(OPTFLAGSC) $(INCDIR) $(DEFS)
CPPCOMP=$(CPP) $(CPPFLAGS) $(OPTFLAGS) $(OPTFLAGSCPP) $(INCDIR) $(DEFS)
FCOMP=$(FC) $(FCFLAGS) $(OPTFLAGS) $(OPTFLAGSF) $(INCDIR) $(DEFS)
F90COMP=$(FC) $(FCFLAGS) $(OPTFLAGS) $(OPTFLAGSF) $(INCDIR) $(DEFS)

ETAGS=$(filter %.o,$(SRC:%.f=$(OBJDIR)/%.o))
ETAGS+=$(filter %.o,$(SRC:%.f90=$(OBJDIR)/%.o))
ETAGS+=$(filter %.o,$(SRC:%.c=$(OBJDIR)/%.o))
ETAGS+=$(filter %.o,$(SRC:%.cpp=$(OBJDIR)/%.o))

# If the BLASLIB is not defined add the included blas to the libs.
# If the LAPACKLIB is not defined add the included lapack to the libs.
FEATLIB+= $(if $(LAPACKLIB),,lapack) $(if $(BLASLIB),,blas)

FEAT=$(FEATLIB:%=$(LIBDIR)/lib%.a)

.PHONY: all
all: greet $(EXEC)
	@echo "Done," $(EXEC) "is ready."
	@echo

.PHONY: greet
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

$(OBJDIR)/%.o : %.cpp
	$(CPPCOMP) -c -o $@ $<

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
	@($(FEATFLOW)/bin/libmkdep.sh $@ $(INCDIR))

.PHONY: clean
clean:
	-rm -f Deps.mk
	-rm -f $(OBJDIR)/*.o 
	-rm -f $(MODDIR)/*.mod 

.PHONY: clean_exec
clean_exec: 
	-rm -f $(EXEC)

.PHONY: purge
purge: clean clean_exec

.PHONY: purge_all
purge_all: purge
	@(for i in *-*-*-* ; do if [ -x $$i ] ; then echo rm $$i ; rm -f $$i ; fi ; done )
	-rm -f -r obj/*
	-rm -f -r \#gmv/* \#avs/* \#points/* \#ns/* \#film/* \#data/\#DX*
	-rm -f -r \#data/*.ERR \#data/*.PRT \#data/*.SYS \#data/*.TRC 

.NOTPARALLEL: clean purge run

debug: all

mpi: all

run: all
	@(time ./$(EXEC))

id: .id
help: .help
-include Deps.mk

##############################################################################
# Documentation
##############################################################################

# Java compiler
JAVAC = javac

# Java runtime engine
JAVA  = java

# Name of documentation parser 
PARSER = p1

##############################################################################
# HTML
##############################################################################

# File names for HTML header and footer
HTML_HEADER_FILE = $(FEATFLOW)/docs/html/header.html
HTML_FOOTER_FILE = $(FEATFLOW)/docs/html/footer.html

HTMLDIR = docs/html

# File name of the resulting HTML documentation
HTML_MASTER_FILE = $(HTMLDIR)/$(shell basename $(shell pwd)).html

.PHONY: html
html: $(HTML_HEADER_FILE) $(HTML_FOOTER_FILE) $(SRC:%=$(HTMLDIR)/%.html)
	@echo "Creating HTML documentation in <$(HTML_MASTER_FILE)>.";
	@rm -f $(HTML_MASTER_FILE)
	@cat $(HTML_HEADER_FILE) >> $(HTML_MASTER_FILE);
        # Creating table of contents
	@echo "<ul>" >> $(HTML_MASTER_FILE);
	@$(foreach file, $(SRC), \
		modulename=$(basename $(file)); \
		echo "<li><a href=\"#module"$${modulename}"\">"$${modulename}"</a></li>" \
		>> $(HTML_MASTER_FILE); )
	@echo "</ul>" >> $(HTML_MASTER_FILE);
	@echo "<hr>" >> $(HTML_MASTER_FILE);
        # Creating documentation
	@$(foreach file, $(SRC:%=$(HTMLDIR)/%.html), \
		modulename=$(basename $(basename $(notdir $(file)))); \
		(echo; echo; \
		echo "<a name=\"module"$${modulename}"\"></a>"; \
		cat ${file}; \
		echo "<div class=\"navigate_to_top\">" ; \
		echo "    <a href=\"#module"$${modulename}"\">To top of this module's documentation</a>"; \
		echo "    <br>"; \
		echo "    <a href=\"#moduleanchor\">To top of documentation</a>"; \
		echo "</div>"; ) >> $(HTML_MASTER_FILE); )
	@cat $(HTML_FOOTER_FILE) >> $(HTML_MASTER_FILE);
	@echo "HTML documentation created and stored in <$(HTML_MASTER_FILE)>.";

.PHONY: ps
ps:
	@echo "Creating documentation (PS) in"
	@pwd
	mkdir -p docs/tex

.PHONY: pdf
pdf:
	@echo "Creating documentation (PDF) in"
	@pwd
	mkdir -p docs/tex

##############################################################################
# Implicit rules
##############################################################################

# Compile Java parser
$(PARSER).class: $(FEATFLOW)/bin/$(PARSER).java
	@echo "Compiling parser...";
	$(JAVAC) -d . $<

# Wrap a Fortran 77 source file to get an XML file
# which we then can pass to the Java-based parser.
%.f.xml: %.f
	@echo; echo "Wrapping $< in XML format...";
	@(echo '<?xml version="1.0" encoding="iso-8859-1" ?>'; \
	  echo '<db>'; \
	  cat $< | sed -f $(FEATFLOW)/bin/prepare4parser.sed; \
	  echo '</db>') > $@;

# Wrap a Fortran 90 source file to get an XML file
# which we then can pass to the Java-based parser.
%.f90.xml: %.f90
	@echo; echo "Wrapping $< in XML format...";
	@(echo '<?xml version="1.0" encoding="iso-8859-1" ?>'; \
	  echo '<db>'; \
	  cat $< | sed -f $(FEATFLOW)/bin/prepare4parser.sed; \
	  echo '</db>') > $@;

# Wrap a C source file to get an XML file
# which we then can pass to the Java-based parser.
%.c.xml: %.c
	@echo; echo "Wrapping $< in XML format...";
	@(echo '<?xml version="1.0" encoding="iso-8859-1" ?>'; \
	  echo '<db>'; \
	  cat $< | sed -f $(FEATFLOW)/bin/prepare4parser.sed; \
	  echo '</db>') > $@;

# Extract module documentation from wrapped FEAT kernel module
# with help of the Java-based parser. 
$(HTMLDIR)/%.html: %.xml $(PARSER).class
	@echo "Parsing $< to create module documentation...";
	@mkdir -p $(HTMLDIR)
	@$(JAVA) -classpath . $(PARSER) html $< $(HTMLDIR);

##############################################################################
# Etags
##############################################################################

.PHONY: tags
tags: $(SRC)
	@(rm -f TAGS)
	@(etags $(filter %.c %.f %.f90,$^))
