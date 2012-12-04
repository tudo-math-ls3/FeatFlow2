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
!# 1.) spdsc_getDist1LevelDiscr
!#     -> Creates a distributed FE discretisation for one spatial level
!#
!# </purpose>
!##############################################################################

module spacediscretisation

  use fsystem
  use storage
  use genoutput
  use basicgeometry
  use cubature
  use element
  use boundary
  use triangulation
  use spatialdiscretisation
  use bcassemblybase
  use collection
  
  use constantsdiscretisation
  use structuresdiscretisation

  implicit none
  
  private
  
  public :: spdsc_getDist1LevelDiscr
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine spdsc_getDist1LevelDiscr (rboundary,rtriangulation,&
      rdiscretisation,rphysics,ieltype)
  
!<description>
  ! Sets up the discretisation structures for the (primal) velocity and pressure
  ! as specified by the parameters in the parameter list rparList.
!</description>

!<input>
  ! Underlying boundary.
  type(t_boundary), intent(IN), target :: rboundary
  
  ! Underlying triangulation
  type(t_triangulation), intent(IN), target :: rtriangulation
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics

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

    select case (rphysics%cequation)
    
    ! *************************************************************
    ! Stokes/Navier Stokes.
    ! *************************************************************
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
    
      ! 2D, 3 equations

      select case (ieltype)
      case (0)
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
        call output_line("Unknown discretisation: ieltype = "//trim(sys_siL(ieltype,10)),&
            OU_CLASS_ERROR,OU_MODE_STD,"spdsc_getDist1LevelDiscr")
        call sys_halt()
      end select
    
    ! *************************************************************
    ! Heat equation
    ! *************************************************************
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
    
      select case (ieltype)
      case (0)
        
        ! Q1
      
        call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                    rtriangulation, rboundary)

        call spdiscr_initDiscr_simple ( &
                    rdiscretisation%RspatialDiscr(1), &
                    EL_Q1,rtriangulation, rboundary)
                    
      case (1)
        
        ! Q2
      
        call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                    rtriangulation, rboundary)

        call spdiscr_initDiscr_simple ( &
                    rdiscretisation%RspatialDiscr(1), &
                    EL_Q2,rtriangulation, rboundary)
                    
      case default
        call output_line("Unknown discretisation: ieltype = "//trim(sys_siL(ieltype,10)),&
            OU_CLASS_ERROR,OU_MODE_STD,"spdsc_getDist1LevelDiscr")
        call sys_halt()
      end select

    case default
      call output_line("Unknown equation",&
          OU_CLASS_ERROR,OU_MODE_STD,"spdsc_getDist1LevelDiscr")
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdsc_getDist1LevelDiscr1D (rboundary,rtriangulation,&
      rdiscretisation,rphysics,ieltype)
  
!<description>
  ! Sets up the discretisation structures for the velocity on the boundary
  ! of the domain.
!</description>

!<input>
  ! Underlying boundary.
  type(t_boundary), intent(IN), target :: rboundary
  
  ! Underlying triangulation
  type(t_triangulation), intent(IN), target :: rtriangulation
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics

  ! Element type (corresponding quad element)
  ! 0 = P1~(E031) / P1~(E031)
  ! 1 = P1~(E030) / P1~(E030)
  ! 2 = P1~(EM31) / P1~(EM31)
  ! 3 = P1~(EM30) / P1~(EM30)
  ! 4 = P2 (E013) / P2 (E013)
  ! 5 = P1~(EM30) / P1~(EM30) unpivoted (much faster than 3 but less stable)
  ! 6 = P1~(EM30) / P1~(EM30) unscaled (slightly faster than 3 but less stable)
  integer, intent(in) :: ieltype
!</input>

!<inputoutput>
  ! The discretisation structure to be set up.
  type(t_blockDiscretisation), intent(inout) :: rdiscretisation
!</inputoutput>

!</subroutine>

!?????????
!    ! local variables
!    integer :: i,nintervals,nel
!    type(t_boundaryRegion) :: rregion
!    
!    ! Determine the number of intervals on the complete boundary.
!    nintervals = 0
!    do i=1,boundary_igetNBoundComp(rboundary)
!      ! Create a region for the complete boundary component
!      call boundary_createRegion (rboundary, i, 0, rregion)
!      
!      ! Number of intervals = number of elements touching that
!      call bcasm_getElementsInBdRegion (rtriangulation,rregion,nel)
!      nintervals = nintervals + nel
!    end do

    select case (rphysics%cequation)
    
    ! *************************************************************
    ! Stokes/Navier Stokes.
    ! *************************************************************
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
    
      ! 2D, 3 equations

      select case (ieltype)
      case (4)
        call spdiscr_initBlockDiscr (rdiscretisation,3,&
                                    rtriangulation, rboundary)
        
        call spdiscr_initDiscr_simple ( &
                    rdiscretisation%RspatialDiscr(1), &
                    EL_Q2,rtriangulation, rboundary)
                    
        call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
            rdiscretisation%RspatialDiscr(2), .true.)
    
        call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
            EL_QP1NP, rdiscretisation%RspatialDiscr(3))
                    
      case default
        call output_line("Unknown discretisation: ieltype = "//trim(sys_siL(ieltype,10)),&
            OU_CLASS_ERROR,OU_MODE_STD,"spdsc_getDist1LevelDiscr1D")
        call sys_halt()
      end select
    
    ! *************************************************************
    ! Heat equation
    ! *************************************************************
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
    
      select case (ieltype)
      case (0)
        
        ! P1
      
        call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                    rtriangulation, rboundary)

        call spdiscr_initDiscr_simple ( &
                    rdiscretisation%RspatialDiscr(1), &
                    EL_Q1,rtriangulation, rboundary)
                    
      case (1)
        
        ! P2
      
        call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                    rtriangulation, rboundary)

        call spdiscr_initDiscr_simple ( &
                    rdiscretisation%RspatialDiscr(1), &
                    EL_Q2,rtriangulation, rboundary)
                    
      case default
        call output_line("Unknown discretisation: ieltype = "//trim(sys_siL(ieltype,10)),&
            OU_CLASS_ERROR,OU_MODE_STD,"spdsc_getDist1LevelDiscr1D")
        call sys_halt()
      end select

    case default
      call output_line("Unknown equation",&
          OU_CLASS_ERROR,OU_MODE_STD,"spdsc_get1LevelDiscrNavSt2D")
      call sys_halt()
    end select

  end subroutine

end module
