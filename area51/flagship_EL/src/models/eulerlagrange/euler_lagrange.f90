!##############################################################################
!# ****************************************************************************
!# <name> euler_lagrange </name>
!# ****************************************************************************


module euler_lagrange

  use afcstabbase
  use bilinearformevaluation
  use boundary
  use boundaryfilter
  use collection
  use derivatives
  use element
  use eulerlagrange_basic
  use eulerlagrange_callback
  use eulerlagrange_callback1d
  use eulerlagrange_callback2d
  use eulerlagrange_callback3d
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use hadaptivity
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocerror
  use pprocgradients
  use pprocindicator
  use pprocsolution
  use problem
  use solveraux
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use thermodynamics
  use timestep
  use timestepaux
  use ucd


   implicit none

   private
   public :: eulerlagrange_init
   public :: eulerlagrange_step
   public :: calculatebarycoords
   public :: findnewelement
   public :: wrongelement
   public :: moveparticle
   public :: checkboundary
   public :: calculatevolumepart



    contains





end module euler_lagrange
