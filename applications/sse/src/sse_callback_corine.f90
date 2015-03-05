!##############################################################################
!# ****************************************************************************
!# <name> sse_callback_corine </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for Corine`s problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in "intf_xxxx.inc" files.
!#
!# </purpose>
!##############################################################################

module sse_callback_corine

  use fsystem
  use storage
  use genoutput
  use derivatives
  use boundary
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use element
  use cubature
  use spatialdiscretisation
  use scalarpde
  use domainintegration
  use collection
  use discretebc
  use discretefbc
  use pprocgradients
  use pprocerror
  
  use sse_base
  use sse_base_corine

  implicit none

  private

end module sse_callback_corine
