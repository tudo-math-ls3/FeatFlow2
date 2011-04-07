!##############################################################################
!# ****************************************************************************
!# <name> cc3d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the cc3d-mini problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) coeff_Stokes
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 1.) coeff_Pressure
!#     -> Returns the coefficients for the pressure matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS_x
!#     -> Returns analytical values for the right hand side of the X-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) coeff_RHS_y
!#     -> Returns analytical values for the right hand side of the Y-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) coeff_RHS_z
!#     -> Returns analytical values for the right hand side of the Z-velocity
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 5.) getBoundaryValues
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

module cc3dmini_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use mprimitives
  use discretebc
  use spatialdiscretisation
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  
  implicit none

contains

! ***************************************************************************

  !<subroutine>

  subroutine coeff_Stokes (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest
    
    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine

! ***************************************************************************
  !<subroutine>
  
  subroutine coeff_Pressure (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest
    
    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! Routine is not called anyway in case of constant coefficients!
    Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

    subroutine coeff_RHS_x (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form corresponding
    ! to the X-velocity.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

    subroutine coeff_RHS_y (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form corresponding
    ! to the Y-velocity.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

    subroutine coeff_RHS_z (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form corresponding
    ! to the Z-velocity.
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!    SUBROUTINE getBoundaryValues (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
!                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
!    
!    USE collection
!    USE spatialdiscretisation
!    USE discretebc
!    
!  !<description>
!    ! This subroutine is called during the discretisation of boundary
!    ! conditions. It calculates a special quantity on the boundary, which is
!    ! then used by the discretisation routines to generate a discrete
!    ! 'snapshot' of the (actually analytic) boundary conditions.
!  !</description>
!    
!  !<input>
!    ! Component specifier.
!    ! For Dirichlet boundary: 
!    !   Icomponents(1) defines the number of the solution component, the value
!    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
!    !   2=2nd solution component, e.g. Y-velocity,...,
!    !   3=3rd solution component, e.g. pressure)
!    ! For pressure drop boundary / normal stress:
!    !   Velocity components that are affected by the normal stress
!    !   (usually "1 2" for x- and y-velocity while returned value musr specify
!    !   the pressure at the boundary)
!    INTEGER, DIMENSION(:), INTENT(in)                           :: Icomponents
!  
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    TYPE(t_spatialDiscretisation), INTENT(in)                   :: rdiscretisation
!    
!    ! Boundary region that is currently being processed.
!    TYPE(t_boundaryRegion), INTENT(in)                          :: rboundaryRegion
!    
!    ! The element number on the boundary which is currently being processed
!    integer, intent(in)                                         :: ielement
!    
!    ! The type of information, the routine should calculate. One of the
!    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
!    ! to return one or multiple information value in the result array.
!    INTEGER, INTENT(in)                                         :: cinfoNeeded
!    
!    ! A reference to a geometric object where information should be computed.
!    ! cinfoNeeded=DISCBC_NEEDFUNC : 
!    !   iwhere = number of the point in the triangulation or
!    !          = 0, if only the parameter value of the point is known; this
!    !               can be found in dwhere,
!    ! cinfoNeeded=DISCBC_NEEDFUNCMID : 
!    !   iwhere = number of the edge in which midpoint the value
!    !            should be computed
!    ! cinfoNeeded=DISCBC_NEEDDERIV : 
!    !   iwhere = number of the point in the triangulation or
!    !          = 0, if only the parameter value of the point is known; this
!    !               can be found in dwhere,
!    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
!    !   iwhere = number of the edge where the value integral mean value
!    !            should be computed
!    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
!    !   iwhere = Number of the edge where the normal stress should be computed.
!    INTEGER(I32), INTENT(in)                                    :: iwhere
!
!    ! A reference to a geometric object where information should be computed.
!    ! cinfoNeeded=DISCBC_NEEDFUNC : 
!    !   dwhere = parameter value of the point where the value should be computed,
!    ! cinfoNeeded=DISCBC_NEEDDERIV : 
!    !   dwhere = parameter value of the point where the value should be computed,
!    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
!    !   dwhere = 0 (not used)
!    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
!    !   dwhere = parameter value of the point on edge iwhere where the normal
!    !            stress should be computed.
!    REAL(DP), INTENT(in)                                        :: dwhere
!     
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    TYPE(t_collection), INTENT(inout), OPTIONAL      :: rcollection
!
!  !</input>
!  
!  !<output>
!    ! This array receives the calculated information. If the caller
!    ! only needs one value, the computed quantity is put into Dvalues(1). 
!    ! If multiple values are needed, they are collected here (e.g. for 
!    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
!    !
!    ! The function may return SYS_INFINITY_DP as a value. This indicates the
!    ! framework to ignore the node and treat it as 'natural boundary condition'
!    ! node.
!    REAL(DP), DIMENSION(:), INTENT(out)                         :: Dvalues
!  !</output>
!    
!  !</subroutine>
!
!    INTEGER :: icomponent
!    REAL(DP) :: x,y
!    
!    REAL(DP), PARAMETER :: dinflowSpeed = 0.3_DP
!    
!    SELECT CASE (cinfoNeeded)
!    CASE (DISCBC_NEEDFUNC,DISCBC_NEEDDERIV,DISCBC_NEEDINTMEAN)
!      ! Get from the current component of the PDE we are discretising:
!      icomponent = Icomponents(1)
!      
!      ! -> 1=X-velocity, 2=Y-velocity.
!      
!      ! Return zero Dirichlet boundary values for all situations by default.
!      Dvalues(1) = 0.0_DP
!
!      ! Now, depending on the problem, calculate the actual velocity value.
!      IF (rboundaryRegion%iboundCompIdx .EQ. 1) THEN
!        SELECT CASE (icomponent)
!        CASE (1) ! X-velocity
!          IF ((dwhere .GE. 3.0_DP) .AND. (dwhere .LE. 4.0_DP)) THEN
!            CALL boundary_getCoords(rdiscretisation%p_rboundary, &
!                          rboundaryRegion%iboundCompIdx, dwhere, &
!                          x, y)
!            Dvalues(1) = mprim_getParabolicProfile (y,0.41_DP,dinflowSpeed)
!          END IF
!
!        CASE (2) ! Y-velocity
!          ! Nothing to do here.
!        END SELECT
!      END IF
!    
!    CASE (DISCBC_NEEDNORMALSTRESS)
!      
!      ! Normal stress on inflow is 0.025, let's say.
!      Dvalues(1) = -0.025_DP
!    
!    END SELECT
!    
!  END SUBROUTINE
!

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues (Icomponents,rdiscretisation,rmeshRegion,&
                                cinfoNeeded,Iwhere,Dwhere,Dcoords,Dvalues,&
                                rcollection)
  
  use collection
  use spatialdiscretisation
  use meshregion
  
!<description>
  ! This subroutine is called during the assembly of boundary conditions which
  ! are defined on mesh regions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the solution component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...,
  !   3=3rd solution component, e.g. pressure)
  ! For pressure drop boundary / normal stress:
  !   Velocity components that are affected by the normal stress
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in)                              :: rmeshRegion

  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
  ! An array holding information about what type of DOF is currently processed.
  ! The information is build up as follows:
  ! Iwhere(1) = vertice number of the DOF, if the DOF is vertice-based, otherwise 0
  ! Iwhere(2) = edge number of the DOF, if the DOF is edge-based, otherwise 0
  ! Iwhere(3) = face number of the DOF, if the DOF is face-based, otherwise 0
  ! Iwhere(4) = currently processed element number.
  ! If Iwhere(1) = Iwhere(2) = Iwhere(3) = 0, then the DOF is element based.
  integer, dimension(4), intent(in)                           :: Iwhere
  
  ! The coordinates of the point which is currently processed, given in
  ! reference coordinates of the currently processed cell type (edge,face,element).
  ! If the DOF is vertice-based, then Dwhere is undefined.
  ! If the DOF is edge-based or element-based in 1D, then Dwhere has dimension 1.
  ! If the DOF is face-based or element-based in 2D, then Dwhere has dimension 2.
  ! IF the DOF is element-based in 3D, then Dwhere has dimension 3.
  real(DP), dimension(:), intent(in)                          :: Dwhere

  ! The coordinates of the point for which the boundary values are to be
  ! calculated.
  real(DP), dimension(:), intent(in)                          :: Dcoords

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  !
  ! The function may return SYS_INFINITY_DP as a value. This indicates the
  ! framework to ignore the node and treat it as 'natural boundary condition'
  ! node.
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

  integer :: icomponent
  real(DP) :: x,y,z
  
  ! Get from the current component of the PDE we are discretising:
  icomponent = Icomponents(1)
  
  ! -> 1=X-velocity, 2=Y-velocity, 3=Z-velocity.
  
    select case (cinfoNeeded)
    case (DISCBC_NEEDFUNC,DISCBC_NEEDDERIV,DISCBC_NEEDINTMEAN)
    
      ! Return zero Dirichlet boundary values for all situations by default.
      Dvalues(1) = 0.0_DP

      ! Now, depending on the problem, calculate the actual velocity value.
      select case (icomponent)
      case (1) ! X-velocity
        x = Dcoords(1)
        y = Dcoords(2)
        z = Dcoords(3)
        if ((x .gt. -0.001_DP) .and. (x .lt. 0.001)) then
          Dvalues(1) = y*(1.0_DP-y)*z*(1.0_DP-z)
        end if

      case (2) ! Y-velocity
        ! Nothing to do here.

      case (3) ! Z-velocity
        ! Nothing to do here.
      end select

    case (DISCBC_NEEDNORMALSTRESS)
      
      ! Normal stress on inflow is 0.025, let's say.
      Dvalues(1) = -0.025_DP
    
    end select
    
  end subroutine

end module
