!##############################################################################
!# ****************************************************************************
!# <name> chemotaxis_callback_blowup </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) coeff_Laplace
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!# 1.2) coeff_chemo
!#     -> Returns the coefficients for the chemo matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
	
!# 2.) coeff_RHS
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getBoundaryValues_c_callback  /  getBoundaryValues_u_callback
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

module chemotaxis_callback_blowup

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use triangulation
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use element
  implicit none
	
    ! Defining pi
    real(DP), parameter, public :: PI = 3.141592654_DP

contains


 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! subroutines for initial prescription of chemo !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  initial_c_callback  !!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine initial_c_callback (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    integer, intent(IN)                                         :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! loop-indices
    integer :: icub, iel

    DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            Dvalues( icub, iel ) = userPresc_chemoattrInitCond( Dpoints(1,icub,iel), &
                                                                Dpoints(2,icub,iel) )
        END DO
    END DO
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! end subroutines for initial prescription of chemo !!!!
  !!!!!!!!!!!!!!!!!!!!!  initial_c_callback  !!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! subroutines for initial prescription of cell !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  initial_u_callback  !!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    subroutine initial_u_callback (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    integer, intent(IN)                                         :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! loop-indices
    integer :: icub, iel

    DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            !2D Dvalues( icub, iel ) =  1.0_DP + ic_pattern  ( Dpoints ( 1, icub, iel ), Dpoints ( 2, icub, iel ) )
            Dvalues( icub, iel ) = userPresc_cellsInitCond( Dpoints(1,icub,iel), & 
                                                            Dpoints(2,icub,iel) )
        END DO
    END DO

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! end subroutines for initial prescription of cell !!!!!
  !!!!!!!!!!!!!!!!!!!!!  initial_u_callback  !!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! subroutines for chemo boundary conditions !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!! getBoundaryValues_c_callback !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine getBoundaryValues_c_callback (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
  !<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
  !</description>
  
  !<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(IN)                                         :: cinfoNeeded
  
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
  integer(I32), intent(IN)                                     :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  real(DP), intent(IN)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(IN), optional      :: rcollection
  !</input>

  !<output>
  ! This array receives the calculated information. 
  real(DP), dimension(:), intent(OUT) :: Dvalues
  
  !</output>

  ! coordinates
  real(DP) :: x,y



  !</subroutine>

  ! assign the corresponding vectorentry of uvector to the needed callback coefficient
  !C_0 = rcollection%DquickAccess(1)
  !U_0 = rcollection%DquickAccess(2)
  
  ! To get the X/Y-coordinates of the boundary point, use:
  CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, x, y)
!     x=Dcoords(1)
!     y=Dcoords(2)

   ! Return zero Dirichlet boundary values for all situations.
   !if ( ((dwhere .ge. 0.0_DP) .and (dwhere .le 1.0_DP)) .or. &
   !     ((dwhere .ge. 2.0_DP) .and (dwhere .le 3.0_DP)) ) rhen
   !    Dvalues(1) = dx*dy
   !end if
   
   Dvalues(1) = userPresc_chemoSol( x, y )

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! end subroutines for chemo boundary conditions !!!!!!!
  !!!!!!!!!!!!!!!!!!!!! getBoundaryValues_c_callback !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! subroutines for cell boundary conditions !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!! getBoundaryValues_u_callback !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine getBoundaryValues_u_callback (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
  !<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
  !</description>
  
  !<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(IN)                                         :: cinfoNeeded
  
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
  integer(I32), intent(IN)                                     :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  real(DP), intent(IN)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(IN), optional      :: rcollection
  !</input>

  !<output>
  ! This array receives the calculated information. 
  real(DP), dimension(:), intent(OUT) :: Dvalues
  
  ! coordinates
  real(DP) :: x,y
!</output>


  !</subroutine>

  ! assign the corresponding vectorentry of uvector to the needed callback coefficient
  !C_0 = rcollection%DquickAccess(1)
  !U_0 = rcollection%DquickAccess(2)
  
  ! To get the X/Y-coordinates of the boundary point, use:
  CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, x, y)
!     x=Dcoords(1)
!     y=Dcoords(2)
   ! Return zero Dirichlet boundary values for all situations.
   !if ( ((dwhere .ge. 0.0_DP) .and (dwhere .le 1.0_DP)) .or. &
   !     ((dwhere .ge. 2.0_DP) .and (dwhere .le 3.0_DP)) ) rhen
   !    Dvalues(1) = dx*dy
   !end if
   
   Dvalues(1) = userPresc_cellsSol( x, y )

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! end subroutines for cell boundary conditions !!!!!!!
  !!!!!!!!!!!!!!!!!!!!! getBoundaryValues_u_callback !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! callback subroutine: rhs for chemoattractant !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  callback_chemo_rfc  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine callback_chemo_rfc (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE    
    integer :: icub, iel
    real(DP) :: x, y
    
    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
!     real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl1, DvaluesFevl2

    ! This is the vector which is of interest
    type(t_vectorScalar):: rvector, ranalytCell

    real(DP) :: PHI
    ! some variables for error-ctrl
    real(DP) :: error_max, error
    
    PHI = rcollection%DquickAccess(1)
    
!     rvector = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
!     ranalytCell = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    ! allocate(DvaluesFevl1(1,npointsPerElement , nelements))
!     allocate(DvaluesFevl2(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
!     call fevl_evaluate_sim4(rvector, &
!                                  rdomainIntSubset, DER_FUNC, DvaluesFevl1, 1)
!     call fevl_evaluate_sim4(ranalytCell, &
!                                  rdomainIntSubset, DER_FUNC, DvaluesFevl2, 1)
    error_max = 0_DP
    ! laplacian
    DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            x=Dpoints(1,icub,iel)
            y=Dpoints(2,icub,iel)
            
            Dcoefficients(1,icub,iel) =  userPresc_chemoSol( x, y ) - &
                                                  PHI*userPresc_cellsSol( x, y )

        END DO
    END DO

     !deallocate(DvaluesFevl1)
!     deallocate(DvaluesFevl2)
 
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! end callback subroutine: rhs for chemoattractant !!!!!!
  !!!!!!!!!!!!!!!!!!!!!  callback_chemo_rfc  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! callback subroutine: rhs for cell density !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  callback_chemo_rfu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine callback_chemo_rfu (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
    real(DP) :: CHI
!</subroutine>

    ! local variables

    ! constants of the PDE
    integer :: icub, iel
    real(DP) :: x, y
    
    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector

    real(DP) :: dtstep

    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)
    
    !coordinates: Dpoints(x,:,:)
        ! laplacian
    DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            x=Dpoints(1,icub,iel)
            y=Dpoints(2,icub,iel)
            ! laplacian
            ! For the linear test-case
            ! eg. u_analytic = (16-x)
            !      c_analytic = x+y+z
!             Dcoefficients(1,icub,iel) = dtstep* (-convecRelaxation)
            ! For the andriy test-case
            ! eg. u_analytic = x*(16-x)
            !      c_analytic = x+y+z
!             Dcoefficients(1,icub,iel) = dtstep*( (2.0_DP +2.0_DP+2.0_DP) / 768.0_DP&
!                                                            + convecRelaxation*(16.0_DP-2.0_DP*x) /18432.0_DP&  
! 							   + convecRelaxation*(16.0_DP-2.0_DP*y) /18432.0_DP&
! 							   + convecRelaxation*(16.0_DP-2.0_DP*z) / 18432.0_DP )

!             Dcoefficients(1,icub,iel) = dtstep*( 2.0_DP &
!                                                            + convecRelaxation*(1.0_DP-2.0_DP*x) )

            Dcoefficients(1,icub,iel) = dtstep*( 2.0_DP &
                                                  - CHI*( ( x-y+2.0_DP*x*y ) / ( x+y )**2.0_DP ) )

        END DO
    END DO
    
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! end callback subroutine: rhs for cell density !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  callback_chemo_rfu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! subroutine for convective term !!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!  callback_K !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine callback_K(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c, rvector_u
    type(t_vectorScalar) :: rvector_x, rvector_y

    ! Some params passed by the collection structure
    real(DP) :: dtstep, CHI, GAMMA, ALPHA, convecRelaxation 

    integer :: icub, iel

    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(3,npointsPerElement,nelements))
    ! Fetching the vector

    rvector_x = collct_getvalue_vecsca (rcollection, "rvector_x",0,'')
    rvector_y = collct_getvalue_vecsca (rcollection, "rvector_y",0,'')
    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)



    ! Fetching the values of rvector_c in the cubature pts.
    call fevl_evaluate_sim4(rvector_x, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_y, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 2)
    !call fevl_evaluate_sim4(rvector_z, &
    !                             rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    !call fevl_evaluate_sim4(rvector_c, &
    !                             rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 1)
    !call fevl_evaluate_sim4(rvector_c, &
    !                             rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 2)
    !call fevl_evaluate_sim4(rvector_c, &
    !                             rdomainIntSubset, DER_DERIV_Z, DvaluesFevl, 3)

    ! These calls are neccessary to fit the signature of f_CHI
    !call fevl_evaluate_sim4(rvector_u, &
    !                             rdomainIntSubset, DER_FUNC, DvaluesFevl, 4)
    call fevl_evaluate_sim4(rvector_c, &
                                rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) = CHI*DvaluesFevl(1,icub,iel)/DvaluesFevl(3,icub,iel)
            ! second term
            Dcoefficients(2,icub,iel) = CHI*DvaluesFevl(2,icub,iel)/DvaluesFevl(3,icub,iel)
            ! third term
            ! Dcoefficients(3,icub,iel) = convecRelaxation*CHI*DvaluesFevl(3,icub,iel)             
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! end subroutine for convective term !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  callback_K !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!! subroutines for construction laplace in the defcorr !
  !!!!!!!!!!!!!!!!!!!!!  callback_defcorr_laplace  !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine callback_defcorr_laplace(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, D_1, N

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(1,npointsPerElement,nelements))

    ! Fetching the vector
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    D_1 = rcollection%DquickAccess(2)
    N = rcollection%DquickAccess(3)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI 

    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the 
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) = dtstep * D_1
            ! second term
            Dcoefficients(2,icub,iel) = dtstep * D_1
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!! end subroutines for construction laplace in the defcorr !!
  !!!!!!!!!!!!!!!!!!!!!  callback_defcorr_laplace  !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!! blow-up restrictive term for the defcorr loop !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  coeff_pattern_growthterm  !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine coeff_pattern_growthterm(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep
    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(1,npointsPerElement,nelements))

    ! Fetching the vector
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI 

    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the 
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            Dcoefficients(1,icub,iel) = DvaluesFevl(1,icub,iel) * ( 1.0_DP - DvaluesFevl(1,icub,iel) )
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!! end blow-up restrictive term for the defcorr loop !!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  coeff_pattern_growthterm  !!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!!!!!!! Subroutine to prescribe the analitycal value of chemo !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!! analyt_c_pattern !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine analyt_c_pattern (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    integer, intent(IN)                                         :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! loop-indices
    integer :: icub, iel

    DO iel = 1, nelements
        DO icub = 1, npointsPerElement             
                !RS: c^2 test-case
!             Dvalues(icub, iel) = Dpoints(1,icub,iel)**2+Dpoints(2,icub,iel)+Dpoints(3,icub,iel)
                !RS: andriy and linear test-case
            Dvalues(icub, iel) = userPresc_chemoattrInitCond( Dpoints(1,icub,iel), & 
                                                            Dpoints(2,icub,iel) )
        END DO
    END DO
  end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!!!!! End of subroutine to prescribe the analitycal value of chemo !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!! analyt_c_pattern !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!! Subroutine to prescribe the analitycal value of cells !!!!!
    !!!!!!!!!!!!!!!!!!!!!!! analyt_u_pattern !!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine analyt_u_pattern (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    integer, intent(IN)                                         :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! loop-indices
    integer :: icub, iel

    DO iel = 1, nelements
        DO icub = 1, npointsPerElement             

            Dvalues( icub, iel ) = userPresc_cellsInitCond( Dpoints(1,icub,iel), & 
                                                            Dpoints(2,icub,iel) )
        END DO
    END DO
  end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!! End of subroutine to prescribe the analitycal value of cells !!!!!
    !!!!!!!!!!!!!!!!!!!!!!! analyt_u_pattern !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function userPresc_cellsSol(x,y) result (func_result)
        ! This function returns the analytical solution for assembling 
        ! the RHSs, BCs ect.
            real(DP) :: x,y
            real(DP) :: func_result

            func_result = 0.0_DP!userPresc_cellsInitCond ( x, y )

        end function userPresc_cellsSol

        function userPresc_chemoSol(x,y) result (func_result)
        ! This function returns the analytical solution for assembling 
        ! the RHSs, BCs ect.
            real(DP) :: x,y
            real(DP) :: func_result

            func_result = 0.0_DP

        end function userPresc_chemoSol
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! User prescribed function for setting initial conditions for cells !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function userPresc_cellsInitCond(x,y) result (func_result)		
	    !
	    ! coordinates
	    real(DP) :: x, y
	    !
	    ! function value
		real(DP) :: func_result

        ! setting initial solution
        func_result = 1000.0_DP * exp ( -100.0_DP* ( ( x - 8.0_DP )**2 + ( y - 8.0_DP)**2)  ) 

	end function userPresc_cellsInitCond


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! User prescribed function for setting initial conditions for cells !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function userPresc_chemoattrInitCond(x,y) result (func_result)		
	    !
	    ! coordinates
	    real(DP) :: x, y
	    !
	    ! function value
		real(DP) :: func_result

        ! setting initial solution
        func_result = 500.0_DP * exp ( -50.0_DP* ( ( x - 8.0_DP )**2 + ( y - 8.0_DP)**2)  ) 

	end function userPresc_chemoattrInitCond 


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!! end subroutines for analytical evalution of chemo !!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_Chemo  !!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine ffunction_Target_Chemo (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax, x, y
    integer :: itimedependence, icub, iel

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            x=Dpoints(1,icub,iel)
            y=Dpoints(2,icub,iel)

            Dvalues(icub,iel) = userPresc_chemoSol( x, y)

        END DO
    END DO
    
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!! end subroutines for analytical evalution of chemo !!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_Chemo  !!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! subroutines for analytical evalution of cell !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_Cells  !!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine ffunction_Target_Cells (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax, x, y, z
    integer :: itimedependence, icub, iel

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            x=Dpoints(1,icub,iel)
            y=Dpoints(2,icub,iel)

            Dvalues(icub,iel) = userPresc_cellsSol( x, y )
        END DO
    END DO

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! subroutines for analytical evalution of cell !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_Cells  !!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! subroutines for analytical evalution of chemo !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_ChemoH1  !!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine ffunction_Target_ChemoH1 (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax, x, y
    integer :: itimedependence, icub, iel

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            x=Dpoints(1,icub,iel)
            y=Dpoints(2,icub,iel)

            IF (cderivative .EQ. DER_FUNC) THEN
                Dvalues(:,:) = userPresc_chemoSol (x, y)
            ELSE IF (cderivative .EQ. DER_DERIV_X) THEN
                Dvalues(:,:) = 1.0_DP
            ELSE IF (cderivative .EQ. DER_DERIV_Y) THEN
                Dvalues(:,:) = 1.0_DP
            END IF

        END DO
    END DO

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! subroutines for analytical evalution of chemo !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_ChemoH1  !!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!! subroutines for analytical evalution of cell !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_CellsH1  !!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine ffunction_Target_CellsH1 (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! If the analytical solution is unknown, this routine does not make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    real(DP) :: dtime,dtimeMax, x, y
    integer :: itimedependence, icub, iel

    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      dtime = rcollection%Dquickaccess(1)
      dtimeMax = rcollection%Dquickaccess(3)
      itimedependence = rcollection%Iquickaccess(1)
    else
      itimedependence = 0
      dtime = 0.0_DP
      dtimeMax = 0.0_DP
    end if

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            x=Dpoints(1,icub,iel)
            y=Dpoints(2,icub,iel)

            IF (cderivative .EQ. DER_FUNC) THEN
                Dvalues(:,:) = userPresc_cellsSol (x, y)
            ELSE IF (cderivative .EQ. DER_DERIV_X) THEN
                Dvalues(:,:) = (1.0_DP - 2*Dpoints(1,icub,iel))
            ELSE IF (cderivative .EQ. DER_DERIV_Y) THEN
                Dvalues(:,:) = (1.0_DP - 2*Dpoints(2,icub,iel))
            END IF

        END DO
    END DO

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!! end subroutines for analytical evalution of cell !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!  ffunction_Target_CellsH1  !!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    

  ! ***************************************************************************
  ! Some subroutine, which I cannnot delete at the moment
  ! ***************************************************************************
  subroutine coeff_hillenX_RHS_c(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: PHI, dtstep
    integer :: icub, iel

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector

    ! Setting the params contained in the collection
    dtstep = rcollection%DquickAccess(1)
    PHI = rcollection%DquickAccess(2)

    ! Fetching the vector    
    rvector = collct_getvalue_vecsca (rcollection, "cbvector",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    
   DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            Dcoefficients(1,icub,iel) = dtstep*PHI*DvaluesFevl(1,icub,iel) 
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine
  ! ***************************************************************************
  ! End of some subroutine, which I cannnot delete at the moment
  ! ***************************************************************************

! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_exp (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    integer, intent(IN)                                         :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  !</output>
    
  !</subroutine>

    ! local variables given by the collection
    real(DP) :: COEFF_A, COEFF_B


    COEFF_A = rcollection%DquickAccess(1)
    COEFF_B = rcollection%DquickAccess(2)
!     B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
! 
! 
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
! 
!       else

    ! The values of an bell-shaped exp fct. are returned
    Dvalues(:,:) =  COEFF_A * dexp ( - COEFF_B *  ( ( Dpoints(1,:,:)  - 0.5_DP )** 2 +  ( Dpoints(2,:,:)  - 0.5_DP ) **2 ) )

  end subroutine

  ! ***************************************************************************

  ! ***************************************************************************
  ! Some subroutine, which I cannnot delete at the moment
  ! ***************************************************************************
  subroutine coeff_cherkurbilf(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration  
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with 
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c, rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, CHI

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector
!     rvector_c => rcollection%p_rvectorQuickAccess1 
!     rvector_u => rcollection%p_rvectorQuickAccess2 
    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI 


    ! Fetching the values of rvector_c in the cubature pts.
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 2)

    ! These calls are neccessary to fit the signature of f_CHI
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 4)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the 
    ! LHS of u_n+1 of the third chertock kurganov example

   DO iel = 1,nelements
       DO icub = 1,npointsPerElement
           ! first term
           Dcoefficients(1,icub,iel) =   CHI * DvaluesFevl(1,icub,iel) 
           ! second term
           Dcoefficients(2,icub,iel) =   CHI * DvaluesFevl(2,icub,iel) 
       END DO
   END DO

    deallocate(DvaluesFevl)

  end subroutine
  ! ***************************************************************************
  ! End of some subroutine, which I cannnot delete at the moment
  ! ***************************************************************************


    ! ***************************************************************************
    ! Function, which I cannnot delete at the moment
    ! ***************************************************************************
    function ic_pattern (x,y) result(f_result)
        implicit none
        intrinsic RANDOM_NUMBER 
        real(DP) :: x, y, f_result, random_num
        if ( sqrt ( (x-8.0_DP)**2 + (y-8.0_DP)**2) <= 2.0_DP ) then 
            CALL RANDOM_NUMBER (random_num)
            !f_result = random_num
            f_result = 1.0_DP + random_num/10.0_DP
            !f_result = rand(0) !1.1 * cos ( 4 * ( PI * sqrt ( (x-8)**2 + (y-8)**2 ) ) / 4 ) **2
        else
            f_result = 1.0_DP
        end if
    end function ic_pattern
    ! ***************************************************************************
    ! End of function, which I cannnot delete at the moment
    ! ***************************************************************************

end module
