#!/usr/bin/env make

########################################################################
#                                                                      #
#                   FINITE ELEMENT ANALYSIS TOOLS 2                    #
#                                                                      #
# Authors: M. Koester, M. Moeller, S. Turek, S. Buijssen               #
#                                                                      #
#                                                                      #
# Contact: Applied Mathematics, TU Dortmund University                 #
#          Vogelpothsweg 87, 44227 Dortmund                            #
#          Germany                                                     #
#                                                                      #
# Web:     http://www.featflow.de/en/software/featflow2.html           #
#          mailto:featflow@featflow.de                                 #
#                                                                      #
########################################################################

# Some variables
SUBDIRS = html tex

# Trunk of file name of resulting documentation
MASTER_FILE = feat2


##############################################################################
# The targets


##############################################################################
# Default target
all:	help

html:
	@$(MAKE) -Chtml $@
	@mv html/$(MASTER_FILE).$@ .

ps:
	@$(MAKE) -Ctex $@
	@mv tex/$(MASTER_FILE).$@ .

pdf:
	@$(MAKE) -Ctex $@
	@mv tex/$(MASTER_FILE).$@ .

clean:
	@$(foreach subdir, $(SUBDIRS), (cd $(subdir) && $(MAKE) $@ ); )

purge:	clean
	@$(foreach subdir, $(SUBDIRS), (cd $(subdir) && $(MAKE) $@ ); )
	-rm -f $(MASTER_FILE).html
	-rm -f $(MASTER_FILE).ps
	-rm -f $(MASTER_FILE).pdf


# print a help screen
help:
	@echo "Usage: make [targets...]"
	@echo
	@echo "where targets include:"
	@echo
	@echo "  help           display this help"
	@echo "  clean          remove working files"
	@echo "  purge          remove documentation and working files"
	@echo
	@echo "  html           Create documentation in HTML format"
	@echo "  ps             Create documentation in Postscript format"
	@echo "  pdf            Create documentation in PDF format"


# For explanation see section 4.6 of GNU Make Manual,
# http://www.gnu.org/software/make/manual/html_mono/make.html#SEC41
.PHONY: all clean html ps pdf
