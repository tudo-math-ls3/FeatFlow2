!##############################################################################
!# ****************************************************************************
!# <Name> zpinch_errorestimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# error estimation for the time-dependent magnetohydrodynamic
!# equations in the one-, two- or three-dimensional domain $\Omega$.
!#
!# The following routines are available:
!#
!# 1.) zpinch_calcAdaptationIndicator
!#     -> Calculates the element-wise indicator for refinement/-coarsening
!#        based on the detectation of shocks and contact discontinuities
!#
!# </purpose>
!##############################################################################

module zpinch_errorestimation

  use fsystem
  use hydro_basic
  use linearsystemblock
  use linearsystemscalar
  use storage
  use triangulation

  implicit none

  private

  public :: zpinch_calcAdaptationIndicator

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcAdaptationIndicator(rsolutionHydro,&
      rsolutionTransport, rindicator)

!<description>
    ! This subroutine computes the element-wise indicator for mesh
    ! refinement and re-coarsening based on the detectation of shocks
    ! and contact discontinuities and the strength of gradient variations
!</description>

!<input>
    ! solution vectors
    type(t_vectorBlock), intent(in) :: rsolutionHydro, rsolutionTransport
!</input>

!<inputoutput>
    ! local feature indicator
    type(t_vectorScalar), intent(inout) :: rindicator
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rdensity, rpressure
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddensity, p_Dpressure, p_Dtracer, p_Dindicator
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisActiveElement
    real(DP) :: dgradient,dlambda_i, dlambda_j
    integer :: iel,jel,ive,nve,i,j,iprotectLayer

    real(DP), parameter :: dEpsS = 0.2_DP
    real(DP), parameter :: dEpsC = 0.1_DP


    ! Set pointer to the underlying triangulation
    p_rtriangulation => rsolutionHydro%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation

    ! Extract primitive variables from conservative solution values
    call hydro_getVariable(rsolutionHydro, 'density', rdensity)
    call hydro_getVariable(rsolutionHydro, 'pressure', rpressure)

    ! Create element-wise indicator
    call lsyssc_createVector(rindicator, p_rtriangulation%NEL, .true.)

    ! Set pointers
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    call lsyssc_getbase_double(rdensity, p_Ddensity)
    call lsyssc_getbase_double(rpressure, p_Dpressure)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lsysbl_getbase_double(rsolutionTransport, p_Dtracer)

    ! Loop over all elements
    elements: do iel = 1, p_rtriangulation%NEL

      ! Get number of vertices at element
      nve = tria_getNVE(p_IverticesAtElement, iel)

      ! Initialize gradient value
      dgradient = 0.0_DP

      ! Check element for shocks and contact discontinuities
      vertices: do ive = 1, nve

        ! Get global vertex numbers
        i = p_IverticesAtElement(ive, iel)
        j = p_IverticesAtElement(mod(ive,nve)+1, iel)

        ! Compute nodal values of tracer
        dlambda_i = p_Dtracer(i)/p_Ddensity(i)
        dlambda_j = p_Dtracer(j)/p_Ddensity(j)

        ! Update maximum gradient function
        dgradient = max(dgradient,&
            abs(abs(p_Ddensity(i))  - abs(p_Ddensity(j)))  /&
              max(abs(p_Ddensity(i)),  abs(p_Ddensity(j))),&
            abs(abs(p_Dpressure(i)) - abs(p_Dpressure(j))) /&
              max(abs(p_Dpressure(i)), abs(p_Dpressure(j))),&
            abs(abs(dlambda_i)      - abs(dlambda_j))      /&
              max(1e-2, abs(dlambda_i), abs(dlambda_j))     )
      end do vertices

      ! If we end up here, then the maximum gradient function is adopted
      p_Dindicator(iel) = dgradient
    end do elements

    ! Release temporal memory
    call lsyssc_releaseVector(rdensity)
    call lsyssc_releaseVector(rpressure)

    ! Add protection layers
    allocate(p_BisActiveElement(p_rtriangulation%NEL))

    do iprotectLayer = 1, 4

      p_BisActiveElement = .false.

      do iel = 1, p_rtriangulation%NEL

        if (p_BisactiveElement(iel)) cycle
        if (p_Dindicator(iel) .le. 0.8) cycle

        do ive = 1, tria_getNVE(p_IverticesAtElement, iel)
          jel = p_IneighboursAtElement(ive, iel)
          if (jel .eq. 0) cycle
          if (p_BisactiveElement(jel)) then
            p_Dindicator(jel) = max(p_Dindicator(jel), p_Dindicator(iel))
          else
            if (p_Dindicator(jel) .lt. 0.8) then
              p_Dindicator(jel) = max(p_Dindicator(jel), p_Dindicator(iel))
              p_BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do
    end do

    deallocate(p_BisActiveElement)

  end subroutine zpinch_calcAdaptationIndicator

end module zpinch_errorestimation
