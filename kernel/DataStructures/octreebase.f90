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
  private

!<constants>

!<constantblock description="Constants for octree structure">

  ! Maximum number of items for each node
  integer, parameter, public :: OTREE_MAX    = 8

  ! Item in "North-West-Front" position
  integer, parameter, public :: OTREE_NWF    = 1

  ! Item in "South-West-Front" position
  integer, parameter, public :: OTREE_SWF    = 2

  ! Item in "South-East-Front" position
  integer, parameter, public :: OTREE_SEF    = 3

  ! Item in "North-East-Front" position
  integer, parameter, public :: OTREE_NEF    = 4

  ! Item in "North-West-Back" position
  integer, parameter, public :: OTREE_NWB    = 5

  ! Item in "South-West-Back" position
  integer, parameter, public :: OTREE_SWB    = 6

  ! Item in "South-East-Back" position
  integer, parameter, public :: OTREE_SEB    = 7

  ! Item in "North-East-Back" position
  integer, parameter, public :: OTREE_NEB    = 8

  ! Position of the status information
  integer, parameter, public :: OTREE_STATUS = 0

  ! Position of the parent information
  integer, parameter, public :: OTREE_PARENT = -1

  ! Position of the "position" of the parent information
  integer, parameter, public :: OTREE_PARPOS = -2

  ! Identifier: Node is empty
  integer, parameter, public :: OTREE_EMPTY  =  0

  ! Identifier: Node is subdivided
  integer, parameter, public :: OTREE_SUBDIV = -1

  ! Identifier: Node is deleted
  integer, parameter, public :: OTREE_DEL    = -2

!</constantblock>

!<constantblock description="Constants for octree bounding-box">

  ! Position of the x-minimal value
  integer, parameter, public :: OTREE_XMIN   =  1

  ! Position of the y-minimal value
  integer, parameter, public :: OTREE_YMIN   =  2

  ! Position of the z-minimal value
  integer, parameter, public :: OTREE_ZMIN   =  3

  ! Position of the x-maximal value
  integer, parameter, public :: OTREE_XMAX   =  4

  ! Position of the y-maximal value
  integer, parameter, public :: OTREE_YMAX   =  5

  ! Position of the z-maximal value
  integer, parameter, public :: OTREE_ZMAX   =  6

!</constantblock>

!<constantblock description="Constants for octree operations">

  ! Operation on octree failed
  integer, parameter, public :: OTREE_FAILED    = -2

  ! Item could not be found in the octree
  integer, parameter, public :: OTREE_NOT_FOUND = -1

  ! Item could be found in the octree
  integer, parameter, public :: OTREE_FOUND     =  0

  ! Item was inserted into the octree
  integer, parameter, public :: OTREE_INSERTED  =  1

  ! Item was deleted from the octree
  integer, parameter, public :: OTREE_DELETED   =  2

!</constantblock>

!</constants>

end module octreebase
