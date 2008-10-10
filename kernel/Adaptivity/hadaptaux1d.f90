!##############################################################################
!# ****************************************************************************
!# <name> hadaptaux1d </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# WARNING: Do not USE this module in your applications unless you really
!#          know what you are doing. This module does no error checking!!!
!#
!# This module contains all auxiliary routines which are required for
!# performing h-adaptivity in 1D. Unlike other modules, all subroutines
!# are declared PUBLIC since they are used by module HADAPTIVITY.
!#
!# The following routines are available:
!#
!#  1.) add_vertex1D
!#      -> Adds a new vertex to the adaptation data structure in 1D
!#
!#  2.) remove_vertex1D
!#      -> Removes an existing vertex from the adaptation data structure in 1D
!#
!#  3.) replace_element1D
!#      -> Replaces an existing element by another element of he same type in 1D
!#
!#  4.) add_element1D
!#      -> Adds a new element to the adaptation data structure in 1D
!#
!#  5.) remove_element1D
!#      -> Removes an existing element from the adaptation data structure in 1D
!#
!#  6.) refine_element1D
!#      -> Refines an element by subdivision into two elements
!#
!#  7.) coarsen_element1D
!#      -> Coarsens two elements by combining them to the macro element
!#
!#  8.) mark_refinement1D
!#      -> Marks elements for refinement in 1D
!#
!#  9.) mark_coarsening1D
!#      -> Marks elements for recoarsening in 1D
!# </purpose>
!##############################################################################

module hadaptaux1d

  use collection
  use fsystem
  use hadaptaux
  
  implicit none
  
  public

end module hadaptaux1d
