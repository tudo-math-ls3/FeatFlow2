#!/usr/bin/env make
########################################################################
#                                                                      #
#         FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T          #
#                                                                      #
# Authors: Ch.Becker,                                                  #
#          S.Buijssen, D.Goeddeke, M.Grajewski, H.Wobker,              #
#          S.Kilian, S.Turek                                           #
#                                                                      #
# Contact: Applied Mathematics, Dortmund University of Technology      #
#          Vogelpothsweg 87, 44227 Dortmund                            #
#          Germany                                                     #
#                                                                      #
# Web:     http://www.feast.tu-dortmund.de                             #
#          mailto:feast@math.tu-dortmund.de                            #
#                                                                      #
########################################################################

# Path to this FEAST installation main directory
ifeq ($(strip $(FEASTBASEDIR)),)
FEASTBASEDIR:=.
endif

include $(FEASTBASEDIR)/Makefile.buildID.inc
include $(FEASTBASEDIR)/Makefile.cpu.inc
