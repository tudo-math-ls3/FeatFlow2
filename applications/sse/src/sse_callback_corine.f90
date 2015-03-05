!##############################################################################
!# ****************************************************************************
!# <name> sse_callback_corine </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for Corine`s problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!#
!#
!# 1.) sse_initMatVec_Corine
!#     -> Initialises the matrices/vectors for the Corine`s problem
!#
!# 2.) sse_initDiscreteBC_Corine
!#     -> Initialises the discrete version of the boundary conditions
!#
!# These callback routines are provided:
!#
!# 1.) coeff_Matrix_Corine
!#     -> Returns analytic valyes for the system matrix.
!#
!# 2.) coeff_RHS_Corine
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 3.) coeff_RHS_Bdr_Corine
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 4.) getBoundaryValues_Corine
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!# </purpose>
!##############################################################################

module sse_callback_corine

  use bcassembly
  use bilinearformevaluation
  use blockmatassemblybase
  use blockmatassemblystdop
  use boundary
  use collection
  use cubature
  use derivatives
  use discretebc
  use discretefbc
  use domainintegration
  use element
  use feevaluation2
  use fsystem
  use genoutput
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use pprocerror
  use pprocgradients
  use scalarpde
  use spatialdiscretisation
  use storage
  use triangulation
  
  use sse_base
  use sse_base_corine

  implicit none

  private

  public :: sse_initMatVec_Corine
  public :: sse_initDiscreteBC_Corine

  public :: coeff_Matrix_Corine
  public :: coeff_RHS_Corine
  public :: coeff_RHS_Bdr_Corine
  public :: getBoundaryValues_Corine
  public :: getReferenceFunction_Corine
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initMatVec_Corine(rproblem)

!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,k

    ! A boundary segment
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    select case(rproblem%cproblemtype)
    case (CORINE_1D)
      
    case (CORINE_2D)
      
    end select
    
  end subroutine sse_initMatVec_Corine

! ***************************************************************************

!<subroutine>

  subroutine sse_initDiscreteBC_Corine(rproblem)

!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_boundaryRegion) :: rboundaryRegion

    select case(rproblem%cproblemtype)
    case (CORINE_1D)
      
    case (CORINE_2D)
      
    end select
  
  end subroutine sse_initDiscreteBC_Corine

! ***************************************************************************

!<subroutine>

  subroutine coeff_Matrix_Corine(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,IdofsTrial,IdofsTest,&
      rdomainIntSubset,Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    integer :: cproblemtype,cproblemsubtype,imatrixpos

    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    imatrixpos      = rcollection%IquickAccess(3)

    select case(cproblemtype)
    case(CORINE_1D) !-----------------------------------------------------------

    case(CORINE_2D) !-----------------------------------------------------------

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_Matrix_Corine")
      call sys_halt()
    end select
    
  end subroutine coeff_Matrix_Corine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Corine(rdiscretisation,rform, &
      nelements,npointsPerElement,Dpoints, &
      IdofsTest,rdomainIntSubset,&
      Dcoefficients,rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

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
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    integer :: cproblemtype,cproblemsubtype,ivectorpos
    
    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    ivectorpos      = rcollection%IquickAccess(3)

    select case(cproblemtype)
    case(CORINE_1D) !-----------------------------------------------------------

    case(CORINE_2D) !-----------------------------------------------------------

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Corine")
      call sys_halt()
    end select

  end subroutine coeff_RHS_Corine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_Corine(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates. According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), allocatable :: Dnx,Dny
    integer :: cproblemtype,cproblemsubtype,ivectorpos,icomp,isegment
    
    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    ivectorpos      = rcollection%IquickAccess(3)
    icomp           = rcollection%IquickAccess(4)
    isegment        = rcollection%IquickAccess(5)

    select case(cproblemtype)
    case(CORINE_1D) !-----------------------------------------------------------

    case(CORINE_2D) !-----------------------------------------------------------

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Corine")
      call sys_halt()
    end select

  end subroutine coeff_RHS_Bdr_Corine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_Corine(Icomponents,rdiscretisation,rboundaryRegion,&
      ielement,cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! "snapshot" of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary:
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(in) :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(in) :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in) :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  integer, intent(in) :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   dwhere = 0 (not used)
  real(DP), intent(in) :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    integer :: cproblemtype,cproblemsubtype,ivectorpos,icomp,isegment
  
    ! Get configuration from quick-access arrays
    cproblemtype    = rcollection%IquickAccess(1)
    cproblemsubtype = rcollection%IquickAccess(2)
    ivectorpos      = rcollection%IquickAccess(3)
    icomp           = rcollection%IquickAccess(4)
    isegment        = rcollection%IquickAccess(5)

    select case(cproblemtype)
    case(CORINE_1D) !-----------------------------------------------------------

    case(CORINE_2D) !-----------------------------------------------------------

    case default !--------------------------------------------------------------
      call output_line("Invalid problem type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"getBoundaryValues_Corine")
      call sys_halt()
    end select

  end subroutine getBoundaryValues_Corine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_Corine(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
      Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in) :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
  
!</subroutine>

  end subroutine getReferenceFunction_Corine
  
end module sse_callback_corine
