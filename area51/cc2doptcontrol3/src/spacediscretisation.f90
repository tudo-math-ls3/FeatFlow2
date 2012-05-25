!##############################################################################
!# ****************************************************************************
!# <name> spacediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic spatial discretisation related routines for
!# CC2D. Here, matrix and RHS creation routines can be found as well as
!# routines to initialise/clean up discretisation structures.
!#
!# The following routines can be found here:
!#
!# 1.) fgetDist1LvDiscrNavSt2D
!#     -> Callback routine for spdsc_getDist1LevelDiscrNavSt2D to be used in
!#        routines from fespacehierarchy.
!#
!# 2.) spdsc_getDist1LevelDiscr
!#     -> Creates a distributed FE discretisation for one spatial level
!#
!# </purpose>
!##############################################################################

module spacediscretisation

  use fsystem
  use storage
  use genoutput
  use cubature
  use element
  use boundary
  use triangulation
  use spatialdiscretisation
  use collection

  implicit none
  
  private
  
  public :: fgetDist1LvDiscrNavSt2D
  public :: spdsc_getDist1LevelDiscr
  
contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine fgetDist1LvDiscrNavSt2D(rtriangulation,rdiscr,rboundary,rcollection)

!<description>
  ! Callback routine wrapper for spdsc_get1LevelDiscrNavSt2D. Allows to create
  ! a discretisation with the routines from fespacehierarchy.
  ! Expects the following information in rcollection:
  !   rcollection%IquickAccess(1) = ieltype
!</description>

  type(t_triangulation), intent(in) :: rtriangulation
  type(t_blockDiscretisation), intent(out) :: rdiscr
  type(t_collection), intent(inout), optional :: rcollection
  type(t_boundary), intent(in), optional :: rboundary
    
!</subroutine>
    
    ! Create the discretisation corresponding to the triangulation
    ! and element type ieltype.
    call spdsc_getDist1LevelDiscr (rboundary,rtriangulation,&
        rdiscr,rcollection%IquickAccess(1))
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdsc_getDist1LevelDiscr (rboundary,rtriangulation,&
      rdiscretisation,ieltype)
  
!<description>
  ! Sets up the discretisation structures for the (primal) velocity and pressure
  ! as specified by the parameters in the parameter list rparList.
!</description>

!<input>
  ! Underlying boundary.
  type(t_boundary), intent(IN), target :: rboundary
  
  ! Underlying triangulation
  type(t_triangulation), intent(IN), target :: rtriangulation
  
  ! Element type.
  ! 0 = Q1~(E031) / Q1~(E031) / Q0
  ! 1 = Q1~(E030) / Q1~(E030) / Q0
  ! 2 = Q1~(EM31) / Q1~(EM31) / Q0
  ! 3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  ! 4 = Q2 (E013) / Q2 (E013) / QP1
  ! 5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  ! 6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  integer, intent(in) :: ieltype
!</input>

!<inputoutput>
  ! The discretisation structure to be set up.
  type(t_blockDiscretisation), intent(inout) :: rdiscretisation
!</inputoutput>

!</subroutine>

    select case (ieltype)
    case (0)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                   rtriangulation, rboundary)

      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_E031,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, rdiscretisation%RspatialDiscr(3))

    case (1)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                  rtriangulation, rboundary)

      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_E030,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, rdiscretisation%RspatialDiscr(3))

    case (2)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                  rtriangulation, rboundary)

      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM31,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, rdiscretisation%RspatialDiscr(3))

    case (3)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                  rtriangulation, rboundary)
      
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, rdiscretisation%RspatialDiscr(3))

    case (4)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                  rtriangulation, rboundary)
      
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_Q2,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_QP1NP, rdiscretisation%RspatialDiscr(3))
                  
    case (5)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                  rtriangulation, rboundary)
      
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30_UNPIVOTED,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, rdiscretisation%RspatialDiscr(3))

    case (6)
      ! Navier Stokes, 2D, 3 equations
      call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                  rtriangulation, rboundary)

      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30_UNSCALED,rtriangulation, rboundary)
                  
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, rdiscretisation%RspatialDiscr(3))

    case default
      call output_line('Unknown discretisation: ieltype = '//trim(sys_siL(ieltype,10)),&
          OU_CLASS_ERROR,OU_MODE_STD,'spdsc_get1LevelDiscrNavSt2D')
      call sys_halt()
    end select
    
  end subroutine

end module
