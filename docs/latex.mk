#!/usr/bin/env make

########################################################################
#                                                                      #
#         FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T          #
#                                                                      #
# Authors: Ch.Becker,                                                  #
#          S.Buijssen, D.Goeddeke, M.Grajewski, H.Wobker,              #
#          S.Kilian, S.Turek                                           #
#                                                                      #
# Contact: Applied Mathematics, University of Dortmund                 #
#          Vogelpothsweg 87, 44227 Dortmund                            #
#          Germany                                                     #
#                                                                      #
# Web:     http://www.feast.uni-dortmund.de                            #
#          mailto:feast@math.uni-dortmund.de                           #
#                                                                      #
########################################################################
#                                                                      #
# Utility Makefile used to create FEAST documentation from LaTeX       #
# sources. Dependencies and external programs are defined below.       #
# To use this file, point the variable TEX_MASTER_FILE in your main    #
# Makefile to the master document (the only one with \begin{document}  #
# in it).                                                              #
#                                                                      #
########################################################################
# $Id: latex.mk,v 2.0 2008/06/01 13:17:50 buijssen Exp $



##############################################################################
# Programs and their flags
#
##############################################################################

# LaTeX and friends
LATEX          = latex
LATEXFLAGS     = 
MAKEINDEX      = makeindex
MAKEINDEXFLAGS =
BIBTEX         = bibtex
BIBTEXFLAGS    =

# Program that generates PostScript files from DVI files.
DVIPS  = dvips
# options when generating PostScript files
DVIPSFLAGS    =
# options when generating PDF files 
DVIPSFLAGSPDF = -Ppdf -G0 -z

# Program that generates PDF files directly from DVI files.
DVIPDF  = dvipdf
# options when generating PDF files 
DVIPDFFLAGS = -sPAPERSIZE=a4 #-Ppdf -G0 -z

# Program that generates PDF files from PostScript files.
PS2PDF = ps2pdf13
# options when generating PDF files 
PS2PDFFLAGS = -sPAPERSIZE=a4

# Program that provides extended grep functionality.
EGREP  = egrep

# Perl
PERL = perl


##############################################################################
# Functions
#
##############################################################################

# Functions to delete a single file / directory at once
# (The reason why we need this and the magic why this works lies within
#  the blank line after the remove command. This way the remove commands
#  in a foreach-loop are triggered one after another - in separate
#  sub-shells. Not in a single sub-shell command. As this may result
#  on some machines in error messages of type:
#  execvp: rm: Arg list too long)
define remove_file
    rm -f $(file)

endef


##############################################################################
# List of dependencies
#
##############################################################################


# All included LaTeX source files
#
# Variant 1 - the easy way: any changed file with extension *.tex triggers
# regeneration of the documentation. #TODO add ../packages and bibtex if */* does 
# not contain these directories already
TEX_DEPS = \
	$(sort $(wildcard *.tex */*.tex ../packages/*.tex ))
BIBTEX_DEPS = \
	$(sort $(wildcard *.bib ../bibtex/*.bib ))
EPS_DEPS = \
	$(sort $(wildcard *.eps */*.eps))
IDX_DEPS = constants functions interfaces modules subroutines types variables


# Variant 2 - a bit more intelligent: parse the master file for \input
# BUG: This has not been adapted to the new directory structure yet
# and \bibliography commands.
#TEX_DEPS = $(shell $(PERL) -ne 'if (s/^[^%]*\s*\\input{([^}]*)}/$$1.tex/) { print; }' $(TEX_MASTER_FILE).tex)
## Bibliographies
#BIBTEX_DEPS = $(shell $(PERL) -ne 'if (s/^[^%]*\s*\\bibliography{([^}]*)}/$$1/) { s/fbmathlit//; s/,/ /g; print; }' $(TEX_MASTER_FILE).tex)
## All included graphics
#EPS_DEPS = \
#	$(sort $(wildcard *.eps */*.eps))


# Variant 3 - more sophisticated way: parse the file for (not commented)
# \input, \includegraphics and \epsfig (deprecated!) commands to find out
# those files that are actually in use. But that requires something like
# a perl script and is considered way too much effort for this purpose.
#
# (not implemented)



##############################################################################
# Implicit rules
#
##############################################################################

# Compile LaTeX document.

%.dvi: MAKEFILE = $(firstword $(MAKEFILE_LIST))
%.dvi %.log %.aux %.toc %.idx: %.tex $(TEX_DEPS) $(BIBTEX_DEPS) $(EPS_DEPS)
	@echo "# $(MAKEFILE): Running $(LATEX) for the first time...";
	$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT);
	@if $(EGREP) "Rerun to get .*references right" $*.log; then \
	    echo "# $(MAKEFILE): Running $(LATEX) again to get references right"; \
	    $(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	fi
	@if $(EGREP) '\\bib(data|cite)' $*.aux; then \
	    echo "# $(MAKEFILE): Making Bibliography using $(BIBTEX)"; \
	    $(BIBTEX) $(BIBTEXFLAGS) $* || exit 1; \
	    if [ -f $*.bbl ]; then \
		echo "# $(MAKEFILE): Running $(LATEX) again to include bibliography"; \
		$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    fi; \
	fi
	makeindex_run=""; \
	$(foreach idxfile, $(IDX_DEPS), \
	    if [ -f $(idxfile)".idx" ]; then \
		echo "# $(MAKEFILE): Creating index using $(MAKEINDEX)"; \
		$(MAKEINDEX) $(MAKEINDEXFLAGS) -o $(idxfile)".ind" $(idxfile)".idx" || exit 1; \
		makeindex_run="1"; \
	    fi;) \
	$(if $$makeindex_run, \
	    echo "# $(MAKEFILE): Running $(LATEX) again to include index"; \
	    $(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT);)
	@-count=5; \
	while $(EGREP) "Rerun to get .*(references|citations) (right|correct)" $*.log && [ $$count -gt 0 ]; do \
	    echo "# $(MAKEFILE): Rerunning $(LATEX), max. $$count times left"; \
	    $(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    count=`expr $$count - 1`; \
	done
	@if [ -f $*.out ] ; then \
	    if $(EGREP) '$*\.out\)' $*.log; then true ; else \
		echo "# $(MAKEFILE): Rerunning $(LATEX) to include PDF outline"; \
		$(LATEX) $(LATEXFLAGS) $< $(LATEXSTDOUT); \
	    fi ; \
	fi
#       print out lists of undefined references/citations
	@if test -n "`$(EGREP) 'Reference .* on page .* undefined on input line .*' $*.log`"; then \
	    echo; \
	    echo "List of undefined references:"; \
	    $(EGREP) 'Reference .* on page .* undefined on input line .*' $*.log; \
	fi
	@if test -n "`$(EGREP) 'Citation .* on page .* undefined on input line .*' $*.log`"; then \
	    echo; \
	    echo "List of undefined citations:"; \
	    $(EGREP) 'Citation .* on page .* undefined on input line .*' $*.log; \
	fi

%.ps : %.dvi
	@if [ -s $< ] ; then \
	    $(DVIPS) $(DVIPSFLAGS) -o $@ $< ; \
	else \
	    echo "# Skipped creating $@: $< does not exist or is empty"; \
	fi

%.pdf : %.dvi
	@if [ ! -s $< ] ; then \
	    echo "# Skipped creating $@: $< does not exist or is empty"; \
	else \
	    path2dvipdf="`which $(DVIPDF)`"; \
	    if [ -x $$path2dvipdf ]; then \
		$(DVIPDF) $(DVIPDFFLAGS) $< $@ ; \
	    else \
		echo "$(DVIPDF) not found, falling back to $(DVIPS) & $(PS2PDF)."; \
		$(DVIPS) $(DVIPSFLAGSPDF) -o- $< | $(PS2PDF) $(PS2PDFFLAGS) - $@ ; \
		rm -f head.tmp body.tmp; \
	    fi; \
	fi

%.clean : %.tex
	-$(foreach file, \
	    $*.aux $*.log $*.toc $*.out \
	    $*.lof $*.lot $*.loa $*.lol $*.thm \
	    $*.nav $*.snm $*.vrb \
	    $*.idx $*.ind $*.ilg $*.glo $*.gls $*.bbl $*.blg, \
	    $(remove_file))
	-$(foreach file, \
	    $(patsubst %, %.idx, $(IDX_DEPS)) \
	    $(patsubst %, %.ilg, $(IDX_DEPS)) \
	    $(patsubst %, %.ind, $(IDX_DEPS)), $(remove_file))

%.purge: %.tex %.clean
	-rm -f $*.pdf $*.ps $*.dvi
        # Emacs backup files
	-rm -f $*.tex~
        # CVS backup files
	-rm -f $*.tex.~[0-9].[0-9]*.~


##############################################################################
# Auxiliary targets
#
##############################################################################

# print a help screen
.PHONY: help
help:
	@echo "Usage: make [targets...]"
	@echo
	@echo "where targets include:"
	@echo
	@echo "  help           display this help"
	@echo
	@echo "  pdf            create documentation in PDF format"
	@echo "  ps             create documentation in Postscript format"
	@echo
	@echo "  clean          remove working files, i.e."
	@echo "                   remove all temporary files LaTeX creates when"
	@echo "                   typesetting a document and which it reads at"
	@echo "                   subsequent invocations. If for some reason,"
	@echo "                   one of these files contains garbage, LaTeX"
	@echo "                   keeps giving error messages. Use this target"
	@echo "                   if you want to start over from scratch"
	@echo "  purge          remove documentation, working and backup files,"
	@echo "                   i.e. in addition to 'make clean' the files"
	@echo "                   $(TEX_MASTER_FILE).pdf, $(TEX_MASTER_FILE).ps, and $(TEX_MASTER_FILE).dvi"
	@echo "                   are removed, along with Emacs and CVS backup"
	@echo "                   files (*~, *.~*~) of the LaTeX sources"


