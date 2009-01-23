!##############################################################################
!# ****************************************************************************
!# <name> afc_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) fcb_coeffRHS
!#     -> compute coefficients for the right-hand side
!#
!# 2.) fcb_setBoundary
!#     -> impose boundary conditions for nonlinear solver
!#
!# </purpose>
!##############################################################################

module afc_callback

  use dofmapping
  use fparser
  use fsystem
  use linearsystemblock
  use linearsystemscalar

  use afc_basic
  use boundaryfilter
  use problem
  use solver

  implicit none

  private
  public :: fcb_coeffRHS
  public :: fcb_coeffTargetFunc
  public :: fcb_setBoundary

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_coeffRHS (rdiscretisation, rform, nelements,&
                           npointsPerElement, Dpoints, IdofsTest,&
                           rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
!</description>
    
!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
!</output>
    
!</subroutine>

    ! local variables
    integer(I32) :: ielement
    
    do ielement = 1, nelements
      call fparser_evalFunction(rrhsParser, 1, 2,&
          Dpoints(:,:,ielement), Dcoefficients(1,:,ielement))
    end do
  end subroutine fcb_coeffRHS

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_coeffTargetFunc (rdiscretisation, rform, nelements,&
                                  npointsPerElement, Dpoints, IdofsTest,&
                                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
!</description>
    
!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
!</output>
    
!</subroutine>

    ! local variables
    integer :: ipoint, ielement
    real(DP) :: x,y

    do ielement = 1, nelements
      do ipoint = 1, npointsPerElement
        x = Dpoints(1, ipoint, ielement)
        y = Dpoints(2, ipoint, ielement)

        if (x .ge. 0.9 .and. x .le. 1.1) then
          Dcoefficients(:, ipoint, ielement) = 1._DP
        else
          Dcoefficients(:, ipoint, ielement) = 0._DP
        end if

!!$        if (sqrt(((x-1.625)**2 + (y-0.25)**2)) .le. 0.125) then
!!$          Dcoefficients(:, ipoint, ielement) = 0.06283185_DP
!!$        else
!!$          Dcoefficients(:, ipoint, ielement) = 0._DP
!!$        end if

!!$        ! Are we in the rectangle [1,2] x [0,0.1]?
!!$        if (x .gt. 1._DP .and. y .le. 0.1) then
!!$          Dcoefficients(:, ipoint, ielement) = 1._DP
!!$        else
!!$          Dcoefficients(:, ipoint, ielement) = 0._DP
!!$        end if
      end do
    end do
    
  end subroutine fcb_coeffTargetFunc
  
  !*****************************************************************************

!<subroutine>

  subroutine fcb_setBoundary(rproblemLevel, rtimestep, rsolver, ru, rres, ru0, istatus)

!<description>
    ! This subroutine imposes the nonlinear boundary conditions.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(IN) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0
!</input>

!<inputoutput>
    ! multigrid level
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: ru

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! OPTIONAL: status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>
    
    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR, &
          NLSOL_PRECOND_NEWTON_FAILED)

      ! Impose nonlinear boundary conditions for solution
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_A),&
                               ru, rres, ru0, rtimestep%dTime)

      
    case (NLSOL_PRECOND_NEWTON)

      ! Impose nonlinear boundary conditions for solution
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
                               ru, rres, ru0, rtimestep%dTime)


    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
                       OU_CLASS_ERROR, OU_MODE_STD,'fcb_setBoundary')
      call sys_halt()
    end select
  end subroutine fcb_setBoundary
  
end module afc_callback
