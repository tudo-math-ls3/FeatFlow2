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

# Path to this FEAST installation main directory
ifeq ($(strip $(FEASTBASEDIR)),)
FEASTBASEDIR:=.
endif

include $(FEASTBASEDIR)/Makefile.buildID.inc
include $(FEASTBASEDIR)/Makefile.cpu.inc
