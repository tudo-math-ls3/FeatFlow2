#!/usr/bin/env make

########################################################################
# FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T  (Release 1.0)   #
#                                                                      #
# Authors: Ch.Becker, S.Kilian, S.Turek                                #
#          Institute of Applied Mathematics & Simulation               #
#          University of Dortmund                                      #
#          D-44227 DORTMUND                                            #
#                                                                      #
# Makefile to create FEAST HTML Documentation from FEAST sources       #
#                                                                      #
# Author: Sven H.M. Buijssen (sven.buijssen@math.uni-dortmund.de)      #
########################################################################

# Name of this Makefile
# (used such that the user can attribute warnings to it)
MAKEFILE = Makefile


##############################################################################
# Paths and invocation of FEAST-specific scripts and settings
#
##############################################################################

# Location of FEAST configure script
CONFIGURE = ../../bin/configure
CONFIGUREFLAGS = --make="$(MAKE) MAKEFLAGS= "

# Location of FEAST's Makefile.inc
MAKEFILE_INC = ../../kernel/arch/Makefile.inc

# Location of FEAST's symbolic link restoration Makefile
MAKEFILE_SYMLINKS = ../../Makefile.symlinks

# List of kernel modules
FEAST_KERNEL_MODULES_WITH_PATH = $(shell $(CONFIGURE) $(CONFIGUREFLAGS) --list-kernel-modules | grep "\.f90" || exit 1;)
FEAST_KERNEL_MODULES           = $(foreach file, $(FEAST_KERNEL_MODULES_WITH_PATH), $(notdir $(file)))


##############################################################################
# Programs and their Flags
#
##############################################################################

# Java compiler
JAVAC = javac

# Java runtime engine
JAVA  = java

# Program to display a line of text
ECHO  = echo


##############################################################################
# Some Variables for creating the documentation from HTML sources
#
##############################################################################

# File name of resulting HTML documentation
HTML_MASTER_FILE = feast.html

# File names for HTML header and footer
HTML_HEADER_FILE = header.html
HTML_FOOTER_FILE = footer.html

# Name of documentation parser 
PARSER = p1

# Suffix for files containing a detailed documentation of a module
SUFFIX_MODULE_DOCUMENTATION = html


##############################################################################
# The targets
#
##############################################################################

# Default target
.PHONY: all
all:	symlinks html

.PHONY: html
html:	| symlinks 
html:	$(HTML_HEADER_FILE) $(HTML_FOOTER_FILE) \
	$(FEAST_KERNEL_MODULES:%=%.$(SUFFIX_MODULE_DOCUMENTATION))
	-rm -f $(HTML_MASTER_FILE);
	@echo "Creating HTML documentation...";
	@cat $(HTML_HEADER_FILE) >> $(HTML_MASTER_FILE);
        # Creating table of contents
	@echo "<ul>" >> $(HTML_MASTER_FILE);
	@$(foreach file, $(FEAST_KERNEL_MODULES), \
		modulename=$(patsubst %.f90,%,$(notdir $(file))); \
		echo "<li><a href=\"#module"$${modulename}"\">"$${modulename}"</a></li>" \
	  >> $(HTML_MASTER_FILE); ) 
	@echo "</ul>" >> $(HTML_MASTER_FILE);
	@echo "<hr>" >> $(HTML_MASTER_FILE);
        # Creating documentation
	@$(foreach file, $(FEAST_KERNEL_MODULES:%=%.$(SUFFIX_MODULE_DOCUMENTATION)), \
		modulename=$(patsubst %.f90.$(SUFFIX_MODULE_DOCUMENTATION),%,$(file)); \
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

# Compile Java parser
$(PARSER).class:	../bin/$(PARSER).java
	@echo "Compiling parser...";
	$(JAVAC) -d . $<

.PHONY: clean
clean:
	-rm -f $(FEAST_KERNEL_MODULES:%=%.$(SUFFIX_MODULE_DOCUMENTATION))
	-rm -f $(PARSER).class $(PARSER)\$$1.class

.PHONY: purge
purge:	clean
	-rm -f $(HTML_MASTER_FILE)


##############################################################################
# Implicit rules
#
##############################################################################

# Wrap a Fortran 90 source file to get an XML file
# which we then can pass to the Java-based parser.
%.f90.xml: %.f90
	@echo; echo "Wrapping $< in XML format...";
	@(echo '<?xml version="1.0" encoding="iso-8859-1" ?>'; \
	  echo '<db>'; \
	  cat $< | sed 's/&//g'; \
	  echo '</db>') > $@;
# Extend search path for Fortran 90 source files to include
# FEAST kernel paths
vpath %.f90 $(sort $(dir $(FEAST_KERNEL_MODULES_WITH_PATH)))

# Extract module documentation from wrapped FEAST kernel module
# with help of the Java-based parser. 
%.$(SUFFIX_MODULE_DOCUMENTATION): %.xml $(PARSER).class
	@echo "Parsing $< to create module documentation...";
	$(JAVA) -classpath . $(PARSER) html $<;


##############################################################################
# Auxiliary targets
#
##############################################################################

.PHONY: symlinks
symlinks:
	@echo "# Creating symbolic links..."
	$(MAKE) -C$(dir $(MAKEFILE_SYMLINKS)) -f $(notdir $(MAKEFILE_SYMLINKS)) symlinks
	@echo "# Symbolic links created."

# print a help screen
.PHONY: help
help:
	@echo "Usage: make [targets...]"
	@echo
	@echo "where targets include:"
	@echo
	@echo "  help           display this help"
	@echo "  clean          remove working files"
	@echo "  purge          remove documentation and working files"
	@echo
	@echo "  html           Create module documentation in HTML format"
	@echo
	@echo "Note:"
	@echo "A FEAST documentation covering more just module descriptions "
	@echo "is only available in PDF or Postscript format."
	@echo "Create it via 'make pdf ps' in directory <../tex>."

# Tell make to delete the target of a rule if it has changed and its commands 
# exit with a nonzero exit status (just as it does when it receives a signal).
# By this rule, make will not create an empty .html file when the java parser
# detects a syntax error in an input file. If we would allow such empty files
# being created, such an error would go undetected in subsequent invocations 
# of make resulting in an incomplete documentation.
.DELETE_ON_ERROR:
