!##############################################################################
!# ****************************************************************************
!# <name> quadtreebase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module file contains basic constants used for quadtrees.
!#
!# </purpose>
!##############################################################################

module quadtreebase

  use fsystem

  implicit none
  private

!<constants>

!<constantblock description="Constants for quadtree structure">

  ! Maximum number of items for each quad
  integer, parameter, public :: QTREE_MAX    = 4

  ! Item in "North-West" position
  integer, parameter, public :: QTREE_NW     = 1

  ! Item in "South-West" position
  integer, parameter, public :: QTREE_SW     = 2

  ! Item in "South-East" position
  integer, parameter, public :: QTREE_SE     = 3

  ! Item in "North-East" position
  integer, parameter, public :: QTREE_NE     = 4

  ! Position of the status information
  integer, parameter, public :: QTREE_STATUS = 0

  ! Position of the parent information
  integer, parameter, public :: QTREE_PARENT = -1

  ! Position of the "position" of the parent information
  integer, parameter, public :: QTREE_PARPOS = -2

  ! Identifier: Quad is empty
  integer, parameter, public :: QTREE_EMPTY  =  0

  ! Identifier: Quad is subdivided
  integer, parameter, public :: QTREE_SUBDIV = -1

  ! Identifier: Quad is deleted
  integer, parameter, public :: QTREE_DEL    = -2

!</constantblock>

!<constantblock description="Constants for quadtree bounding-box">

  ! Position of the x-minimal value
  integer, parameter, public :: QTREE_XMIN   =  1

  ! Position of the y-minimal value
  integer, parameter, public :: QTREE_YMIN   =  2

  ! Position of the x-maximal value
  integer, parameter, public :: QTREE_XMAX   =  3

  ! Position of the y-maximal value
  integer, parameter, public :: QTREE_YMAX   =  4

!</constantblock>

!<constantblock description="Constants for quadtree operations">

  ! Operation on quadtree failed
  integer, parameter, public :: QTREE_FAILED       = -2

  ! Item could not be found in the quadtree
  integer, parameter, public :: QTREE_NOT_FOUND    = -1

  ! Item could be found in the quadtree
  integer, parameter, public :: QTREE_FOUND        =  0

  ! Item was inserted into the quadtree
  integer, parameter, public :: QTREE_INSERTED     =  1

  ! Item was deleted from the quadtree
  integer, parameter, public :: QTREE_DELETED      =  2

  ! Item was moved in the quadtree
  integer, parameter, public :: QTREE_REPOSITIONED =  3

!</constantblock>

!</constants>

end module quadtreebase
