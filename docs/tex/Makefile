#!/usr/bin/env make

########################################################################
# FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T  (Release 1.0)   #
#                                                                      #
# Authors: Ch.Becker, S.Kilian, S.Turek                                #
#          Institute of Applied Mathematics & Simulation               #
#          University of Dortmund                                      #
#          D-44227 DORTMUND                                            #
#                                                                      #
# Makefile to create FEAT2 Documentation from FEAT2 and LaTeX sources  #
#                                                                      #
# Author: Sven H.M. Buijssen (sven.buijssen@math.uni-dortmund.de)      #
########################################################################

# Name of this Makefile
# (used such that the user can attribute warnings to it)
MAKEFILE = Makefile

# Basic FEATFLOW directory;
# define the root directory of the Featflow installation
FEATFLOW=../..


##############################################################################
# Paths and invocation of FEAT2-specific scripts and settings
#
##############################################################################

include ../../kernel/kernel.mk

# Location of FEAST's Makefile.inc
MAKEFILE_INC = ../../kernel/arch/Makefile.inc

# Location of FEAST's symbolic link restoration Makefile
MAKEFILE_SYMLINKS = ../../Makefile.symlinks

# List of kernel modules
#KERNEL_MODULES_WITH_PATH = $(shell $(CONFIGURE) $(CONFIGUREFLAGS) --list-kernel-modules | grep "\.f90" || exit 1;)
KERNEL_MODULES           = $(KERNELSRC)


##############################################################################
# Programs and their Flags
#
##############################################################################

# Java compiler
JAVAC = javac

# Java runtime engine
JAVA = java

# LaTeX
LATEX = latex
#LATEXFLAGS= -interaction=\nonstopmode # -file-line-error-style 

# PDFLaTeX
PDFLATEX      = pdflatex
PDFLATEXFLAGS = $(LATEXFLAGS)

# Program used for index creation.
MAKEINDEX = makeindex
#MAKEINDEXFLAGS ?= -c -g -s gind.ist

# Specifies the bibliography generation program.
BIBTEX = bibtex

# Program that generates PostScript files from DVI files.
DVIPS  = dvips
# options when generating PostScript files
DVIPSFLAGS    =
# options when generating PDF files 
DVIPSFLAGSPDF = -Ppdf -G0

# Program that generates PDF files from PostScript files.
PS2PDF = ps2pdf13

# Program that provides extended grep functionality.
EGREP  = egrep

# Program to display a line of text
ECHO   = echo

##############################################################################
# Some Variables for creating the documentation from LaTeX sources
#
##############################################################################

# Name of master file for documentation
TEX_MASTER_FILE = feat2

# Main dependencies of the LaTeX file
INCLUDES = $(shell $(EGREP) "\include{" $(TEX_MASTER_FILE).tex)
DEPS     = $(INCLUDES:\include{%}=%.tex)

# Name of documentation parser 
PARSER = p1

# Name of LaTeX file that includes short descriptions of each module
TEX_MASTER_FILE_MODULE_PURPOSE = modules_overview_dynamic.tex

# Name of LaTeX file that includes detailed documentation of each module
TEX_MASTER_FILE_MODULE_DOCUMENTATION = modules_description_dynamic.tex

# Suffix for files containing a short description of a module
SUFFIX_MODULE_PURPOSE = descr.tex

# Suffix for files containing a detailed documentation of a module
SUFFIX_MODULE_DOCUMENTATION = docu.tex


##############################################################################
# The targets
#
##############################################################################

# Default target
.PHONY: all
all:	symlinks pdf

.PHONY: tex
tex:	| symlinks
tex:	$(TEX_MASTER_FILE_MODULE_PURPOSE) \
	$(TEX_MASTER_FILE_MODULE_DOCUMENTATION)

dvi:	tex $(TEX_MASTER_FILE).dvi

.PHONY: pdf
pdf:	dvi $(TEX_MASTER_FILE).pdf

.PHONY: ps
ps:	dvi $(TEX_MASTER_FILE).ps

# Create a LaTeX file that includes the descriptions
# for all FEAST kernel modules (via LaTeX \include{}).
$(TEX_MASTER_FILE_MODULE_PURPOSE):	$(KERNEL_MODULES:%=%.$(SUFFIX_MODULE_PURPOSE))
	-rm -f $(TEX_MASTER_FILE_MODULE_PURPOSE);
	@echo "Creating LaTeX file with references to module descriptions...";
	@$(foreach file, $(KERNEL_MODULES), \
		texfile=$(addsuffix .$(SUFFIX_MODULE_PURPOSE), $(file)); \
		modulename=$(patsubst %.f90,%,$(notdir $(file))); \
		(echo "\subsubsection{Module "$${modulename}"}"; \
		 cat $${texfile}; \
		 echo) >> $(TEX_MASTER_FILE_MODULE_PURPOSE); )

# Create a LaTeX file that includes detailed documentation
# for all FEAST kernel modules (via LaTeX \include{}).
$(TEX_MASTER_FILE_MODULE_DOCUMENTATION):	$(KERNEL_MODULES:%=%.$(SUFFIX_MODULE_DOCUMENTATION))
	-rm -f $(TEX_MASTER_FILE_MODULE_DOCUMENTATION);
	@echo "Creating LaTeX file with references to module documentation...";
	@$(foreach file, $(KERNEL_MODULES), \
		texfile=$(addsuffix .$(SUFFIX_MODULE_DOCUMENTATION), $(file)); \
		modulename=$(patsubst %.f90,%,$(notdir $(file))); \
		(cat $${texfile}; \
		 echo) >> $(TEX_MASTER_FILE_MODULE_DOCUMENTATION); )

# Compile Java parser
$(PARSER).class:	../bin/$(PARSER).java
	@echo "Compiling parser...";
	$(JAVAC) -d . $<

.PHONY: clean
clean:	$(TEX_MASTER_FILE).clean
	-rm -f $(DEPS:%.tex=%.aux)
	-rm -f $(addsuffix .xml, $(KERNEL_MODULES));
	-rm -f $(addsuffix .$(SUFFIX_MODULE_PURPOSE), $(KERNEL_MODULES));
	-rm -f $(addsuffix .$(SUFFIX_MODULE_DOCUMENTATION), $(KERNEL_MODULES));
	-rm -f $(TEX_MASTER_FILE_MODULE_PURPOSE);
	-rm -f $(TEX_MASTER_FILE_MODULE_DOCUMENTATION);
	-rm -f $(PARSER).class $(PARSER)\$$1.class

.PHONY: purge
purge:	$(TEX_MASTER_FILE).clean clean $(TEX_MASTER_FILE).purge


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
	  cat $< | sed -e 's/&//g' -e 's/<=/\&lt;=/g' -e 's/<-/\&lt;-/g' -e 's/<[^a-zA-Z?!\/]/\&lt; /g'; \
	  echo '</db>') > $@;
# Extend search path for Fortran 90 source files to include
# FEAST kernel paths
#vpath %.f90 $(sort $(dir $(KERNEL_MODULES_WITH_PATH)))

# Extract module description from wrapped FEAST kernel module
# with help of the Java-based parser. 

%.$(SUFFIX_MODULE_PURPOSE): %.xml $(PARSER).class
	@echo; echo "Parsing $< for module description...";
	$(JAVA) -classpath . $(PARSER) texmod $<;

# Extract module documentation from wrapped FEAST kernel module
# with help of the Java-based parser. 
%.$(SUFFIX_MODULE_DOCUMENTATION): %.xml $(PARSER).class
	@echo; echo "Parsing $< to create module documentation...";
	$(JAVA) -classpath . $(PARSER) tex $<;

# Compile LaTeX document.
%.dvi %.log %.aux %.toc %.idx: %.tex $(DEPS)
	@$(ECHO) "$(MAKEFILE): Running $(LATEX) for the first time...";
	$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT);
	@if $(EGREP) "Rerun to get .*references right" $*.log; then \
	    $(ECHO) "$(MAKEFILE): Running $(LATEX) again to get references right"; \
	    $(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	fi
	@if $(EGREP) '\\bib(data|cite)' $*.aux; then \
	    $(ECHO) "$(MAKEFILE): Making Bibliography using $(BIBTEX)"; \
	    $(BIBTEX) $(BIBTEXFLAGS) $*; \
	    if [ -f $*.bbl ]; then \
		$(ECHO) "$(MAKEFILE): Running $(LATEX) again to include bibliography"; \
		$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    fi; \
	fi
	if [ -f $*.idx ]; then \
	    $(ECHO) "$(MAKEFILE): Creating index using $(MAKEINDEX)"; \
	    $(MAKEINDEX) $(MAKEINDEXFLAGS) -o $*.ind $*.idx; \
	    if [ -s $*.ind ]; then \
		$(ECHO) "$(MAKEFILE): Running $(LATEX) again to include index"; \
		$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    fi; \
	fi
	@-count=5; \
	while $(EGREP) "Rerun to get .*(references|citations) (right|correct)" $*.log && [ $$count -gt 0 ]; do \
	    $(ECHO) "$(MAKEFILE): Rerunning $(LATEX), max. $$count times left"; \
	    $(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    count=`expr $$count - 1`; \
	done
	if [ -f $*.out ] ; then \
	    if $(EGREP) '$*\.out\)' $*.log; then true ; else \
		$(ECHO) "$(MAKEFILE): Rerunning $(LATEX) to include PDF outline"; \
	    	$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    fi ; \
	fi

%.ps : %.dvi
	@if [ -s $< ] ; then \
	    $(DVIPS) $(DVIPSFLAGS) -o $@ $< ; \
	else \
	    $(ECHO) Skipped creating $@: $< does not exist or is empty; \
	fi

%.pdf : %.dvi
	if [ -s $< ] ; then \
	    $(DVIPS) $(DVIPSFLAGSPDF) -o- $< | $(PS2PDF) $(PS2PDFFLAGS) - $@ ; \
	else \
	    $(ECHO) Skipped creating $@: $< does not exist or is empty; \
	fi

%.clean : %.tex
	-rm -f $*.aux $*.log $*.toc $*.out
	-rm -f $*.lof $*.lot $*.loa $*.lol $*.thm
	-rm -f $*.nav $*.snm $*.vrb
	-rm -f $*.idx $*.ind $*.ilg $*.glo $*.gls $*.bbl $*.blg

%.purge: %.tex %.clean
	-rm -f $*.pdf $*.ps $*.dvi


##############################################################################
# Auxiliary targets
#
##############################################################################

.PHONY: 	symlinks
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
#	@echo "  dvi            Create documentation in DVI format"
	@echo "  ps             Create documentation in Postscript format"
	@echo "  pdf            Create documentation in PDF format"

# Tell make to delete the target of a rule if it has changed and its commands 
# exit with a nonzero exit status (just as it does when it receives a signal).
# By this rule, make will not create an empty .tex file when the java parser
# detects a syntax error in an input file. If we would allow such empty files
# being created, such an error would go undetected in subsequent invocations 
# of make resulting in an incomplete documentation.
.DELETE_ON_ERROR:
