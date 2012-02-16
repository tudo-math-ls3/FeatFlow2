!##############################################################################
!# ****************************************************************************
!# <name> octreebase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module file contains basic constants used for octrees.
!#
!# </purpose>
!##############################################################################

module octreebase

  use fsystem

  implicit none

!<constants>

!<constantblock description="Constants for octree structure">
  
  ! Maximum number of items for each node
  integer, parameter :: OTREE_MAX    = 8

  ! Item in "North-West-Front" position
  integer, parameter :: OTREE_NWF    = 1

  ! Item in "South-West-Front" position
  integer, parameter :: OTREE_SWF    = 2

  ! Item in "South-East-Front" position
  integer, parameter :: OTREE_SEF    = 3

  ! Item in "North-East-Front" position
  integer, parameter :: OTREE_NEF    = 4

  ! Item in "North-West-Back" position
  integer, parameter :: OTREE_NWB    = 5

  ! Item in "South-West-Back" position
  integer, parameter :: OTREE_SWB    = 6

  ! Item in "South-East-Back" position
  integer, parameter :: OTREE_SEB    = 7

  ! Item in "North-East-Back" position
  integer, parameter :: OTREE_NEB    = 8
  
  ! Position of the status information
  integer, parameter :: OTREE_STATUS = 0
  
  ! Position of the parent information
  integer, parameter :: OTREE_PARENT = -1

  ! Position of the "position" of the parent information
  integer, parameter :: OTREE_PARPOS = -2

  ! Identifier: Node is empty
  integer, parameter :: OTREE_EMPTY  =  0

  ! Identifier: Node is subdivided
  integer, parameter :: OTREE_SUBDIV = -1

  ! Identifier: Node is deleted
  integer, parameter :: OTREE_DEL    = -2

!</constantblock>

!<constantblock description="Constants for octree bounding-box">

  ! Position of the x-minimal value
  integer, parameter :: OTREE_XMIN   =  1

  ! Position of the y-minimal value
  integer, parameter :: OTREE_YMIN   =  2

  ! Position of the z-minimal value
  integer, parameter :: OTREE_ZMIN   =  3

  ! Position of the x-maximal value
  integer, parameter :: OTREE_XMAX   =  4

  ! Position of the y-maximal value
  integer, parameter :: OTREE_YMAX   =  5

  ! Position of the z-maximal value
  integer, parameter :: OTREE_ZMAX   =  6

!</constantblock>
  
!<constantblock description="Constants for octree operations">

  ! Operation on octree failed
  integer, parameter :: OTREE_FAILED    = -2

  ! Item could not be found in the octree
  integer, parameter :: OTREE_NOT_FOUND = -1

  ! Item could be found in the octree
  integer, parameter :: OTREE_FOUND     =  0

  ! Item was inserted into the octree
  integer, parameter :: OTREE_INSERTED  =  1

  ! Item was deleted from the octree
  integer, parameter :: OTREE_DELETED   =  2
  
!</constantblock>

!</constants>

end module octreebase
