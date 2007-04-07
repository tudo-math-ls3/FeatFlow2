!#########################################################################
!# ***********************************************************************
!# <name> pprocerrror </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating errors of
!# finite element functions.
!#
!# The following routines can be found in this module:
!#
!# 1.) pperr_scalar
!#     -> Calculate $L_2$-error or $H_1$ error to an analytic reference
!#        function.
!# </purpose>
!#########################################################################

MODULE pprocerror

  USE fsystem
  USE storage
  USE boundary
  USE cubature
  USE triangulation
  USE linearsystemscalar
  USE linearsystemblock
  USE spatialdiscretisation
  USE domainintegration
  USE collection
  USE feevaluation

  IMPLICIT NONE

!<constants>

!<constantblock description = "Identifiers for the type of error to be computed.">

  ! $L_2$-error
  INTEGER, PARAMETER :: PPERR_L2ERROR = 1
  
  ! $H_1$-error
  INTEGER, PARAMETER :: PPERR_H1ERROR = 2
  
!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! Number of elements to handle simultaneously when building vectors
  INTEGER :: PPERR_NELEMSIM   = 1000
  
!</constantblock>

!</constants>

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalar (rvectorScalar,cerrortype,derror,&
                           ffunctionReference,rcollection,rdiscretisation)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the cubature
  ! formula to use for each element distribution.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN), TARGET :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! A callback function that provides the analytical reference function
  ! to which the error shouold be computed.
  INCLUDE 'intf_refFunctionSc.inc'
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  TYPE(t_collection), INTENT(IN), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! If not specified, the discretisation structure in the vector is used.
  ! If specified, the discretisation structure must be 'compatible' to the
  ! vector (concerning NEQ,...). pperr_scalar uses the cubature formula
  ! specifier of the linear form in rdiscretisation to compute the integrals
  ! for the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET, OPTIONAL :: rdiscretisation
!</input>

!<output>
  ! Array receiving the calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    TYPE(t_collection), POINTER :: p_rcollection
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! Let p_rcollection point to rcollection - or NULL if it's not
    ! given.
    IF (PRESENT(rcollection)) THEN
      p_rcollection => rcollection
    ELSE
      p_rcollection => NULL()
    END IF
    
    ! Get the correct discretisation structure and check if we can use it.
    IF (PRESENT(rdiscretisation)) THEN
      p_rdiscretisation => rdiscretisation
      CALL lsyssc_checkDiscretisation (rvectorScalar,p_rdiscretisation)
    ELSE
      p_rdiscretisation => rvectorScalar%p_rspatialdiscretisation
    END IF
    
    IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
      PRINT *,'pperr_scalar: No discretisation structure!'
      STOP
    END IF
    
    IF (p_rdiscretisation%ndimension .NE. NDIM2D) THEN
      PRINT *,'pperr_scalar: 3D discretisation currently not supported.'
      STOP
    END IF

    ! The vector must be unsorted, otherwise we can't set up the vector.
    IF (rvectorScalar%isortStrategy .GT. 0) THEN
      PRINT *,'pperr_scalar: Vector must be unsorted!'
      STOP
    END IF
  
    ! Do we have a uniform triangulation? Would simplify a lot...
    IF (p_rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN 
    
      IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
        CALL pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,&
                           ffunctionReference,p_rcollection,p_rdiscretisation)
      ELSE
        PRINT *,'pperr_scalar: Single precision vectors currently not supported!'
      END IF
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    ELSE IF (p_rdiscretisation%ccomplexity .EQ. SPDISC_CONFORMAL) THEN 
    
      IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
        CALL pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,&
                           ffunctionReference,p_rcollection,p_rdiscretisation)
      ELSE
        PRINT *,'pperr_scalar: Single precision vectors currently not supported!'
      END IF
    
    ELSE
      PRINT *,'pperr_scalar: General discretisation not implemented!'
      STOP
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,ffunctionReference,&
                                  p_rcollection,rdiscretisation)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN), TARGET :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! A callback function that provides the analytical reference function
  ! to which the error shouold be computed.
  INCLUDE 'intf_refFunctionSc.inc'

  ! A discretisation structure specifying how to compute the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! A pointer to the collection structure to pass to the callback routine;
  ! or NULL, if none such a structure exists.
  TYPE(t_collection), POINTER            :: p_rcollection
!</input>

!<output>
  ! Array receiving the calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: i,j,k,icurrentElementDistr, ICUBP, NVE
    LOGICAL :: bnonparTrial
    INTEGER(I32) :: IEL, IELmax, IELset
    REAL(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
    
    ! Cubature point coordinates on the reference element
    REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
    
    ! number of cubature points on the reference element
    INTEGER :: ncubp
    
    ! Number of local degees of freedom for test functions
    INTEGER :: indofTrial
    
    ! The triangulation structure - to shorten some things...
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! A pointer to an element-number list
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

    ! An array receiving the coordinates of cubature points on
    ! the real element for all elements in a set.
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal

    ! Pointer to the point coordinates to pass to the element function.
    ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
    ! the trial element is parametric or not.
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
    
    ! Array with coordinates of the corners that form the real element.
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
    
    ! Arrays for saving Jacobian determinants and matrices
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
    
    ! Pointer to KVERT of the triangulation
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
    ! Pointer to DCORVG of the triangulation
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
    
    ! Current element distribution
    TYPE(t_elementDistribution), POINTER :: p_elementDistribution
    
    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    TYPE(t_domainIntSubset) :: rintSubset
    
    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial
  
    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .FALSE.
    SELECT CASE (cerrortype)
    CASE (PPERR_L2ERROR) 
      Bder(DER_FUNC) = .TRUE.
    CASE (PPERR_H1ERROR) 
      Bder(DER_DERIV_X) = .TRUE.
      Bder(DER_DERIV_Y) = .TRUE.
    CASE DEFAULT
      PRINT *,'pperr_scalar2d_conf: Unknown error type identifier!'
      STOP
    END SELECT
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocaing some arrays.
    nelementsPerBlock = MIN(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Get a pointer to the KVERT and DCORVG array
    CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                               p_IverticesAtElement)
    CALL storage_getbase_double2D(p_rtriangulation%h_DcornerCoordinates, &
                               p_DcornerCoordinates)

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_elementDistribution => rdiscretisation%RelementDistribution(icurrentElementDistr)
    
      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_elementDistribution%itestElement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_elementDistribution%itrialElement)
      
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_elementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      j = elem_igetCoordSystem(p_elementDistribution%itrialElement)
      
      ! Allocate memory and get local references to it.
      CALL domint_initIntegration (rintSubset,nelementsPerBlock,ncubp,j,&
          p_rtriangulation%ndim,NVE)
      p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
      p_DcubPtsReal => rintSubset%p_DcubPtsReal
      p_Djac =>        rintSubset%p_Djac
      p_Ddetj =>       rintSubset%p_Ddetj
      p_Dcoords =>     rintSubset%p_DCoords

      ! Put the cubature point coordinates in the right format to the
      ! cubature-point array.
      ! Initialise all entries in p_DcubPtsRef with the same coordinates -
      ! as the cubature point coordinates are identical on all elements
      DO j=1,SIZE(p_DcubPtsRef,3)
        DO i=1,ncubp
          DO k=1,SIZE(p_DcubPtsRef,1)
            ! Could be solved using the TRANSPOSE operator - but often is's 
            ! faster this way...
            p_DcubPtsRef(k,i,j) = Dxi(i,k)
          END DO
        END DO
      END DO
      
      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,4))
    
      ! Check if one of the trial/test elements is nonparametric
      bnonparTrial  = elem_isNonparametric(p_elementDistribution%itestElement)
                      
      ! Let p_DcubPtsTest point either to p_DcubPtsReal or
      ! p_DcubPtsRef - depending on whether the space is parametric or not.
      IF (bnonparTrial) THEN
        p_DcubPtsTrial => p_DcubPtsReal
      ELSE
        p_DcubPtsTrial => p_DcubPtsRef
      END IF
      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Set the current error to 0 and add the error contributions of each element
      ! to that.
      derror = 0.0_DP
                                
      ! Loop over the elements - blockwise.
      DO IELset = 1, p_rtriangulation%NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(p_rtriangulation%NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     .TRUE.,IdofsTrial)
                                     
        ! We have the coordinates of the cubature points saved in the
        ! coordinate array from above. Unfortunately for nonparametric
        ! elements, we need the real coordinate.
        ! Furthermore, we anyway need the coordinates of the element
        ! corners and the Jacobian determinants corresponding to
        ! all the points.
        !
        ! At first, get the coordinates of the corners of all the
        ! elements in the current set. 
        
!        DO IEL=1,IELmax-IELset+1
!          p_Dcoords(:,:,IEL) = p_DcornerCoordinates(:, &
!                              p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!        END DO
!        DO IEL=1,IELmax-IELset+1
!          DO J = 1,NVE
!            DO I = 1,NDIM2D
!              p_Dcoords(I,J,IEL) = p_DcornerCoordinates(I, &
!                                 p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
!            END DO
!          END DO
!        END DO
        CALL trafo_getCoords_sim (elem_igetTrafoType(p_elementDistribution%itrialElement),&
            p_rtriangulation,p_IelementList(IELset:IELmax),p_Dcoords)
        
        ! Depending on the type of transformation, we must now choose
        ! the mapping between the reference and the real element.
        ! In case we use a nonparametric element as test function, we need the 
        ! coordinates of the points on the real element, too.
        ! Unfortunately, we need the real coordinates of the cubature points
        ! anyway for the function - so calculate them all.
        CALL trafo_calctrafo_sim (&
              rdiscretisation%RelementDistribution(icurrentElementDistr)%ctrafoType,&
              IELmax-IELset+1,ncubp,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1),&
              p_DcubPtsReal)
      
        ! Prepare the call to the evaluation routine of the analytic function.    
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
    
        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        SELECT CASE (cerrortype)
        
        CASE (PPERR_L2ERROR)
        
          ! L2-error uses only the values of the function.
          !
          ! It's time to call our coefficient function to calculate the
          ! function values in the cubature points:  u(x,y)
          ! The result is saved in Dcoefficients(:,:,1)
          
          CALL ffunctionReference (DER_FUNC,rdiscretisation, &
                      INT(IELmax-IELset+1),ncubp,p_DcubPtsReal, &
                      IdofsTrial,rintSubset,p_rcollection, &
                      Dcoefficients(:,1:IELmax-IELset+1_I32,1))
          
          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
               p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
               p_elementDistribution%itrialElement, IdofsTrial, &
               ncubp, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_FUNC,&
               Dcoefficients(:,1:IELmax-IELset+1_I32,2))        
          
          ! Subtracttion of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u-u_h,u-u_h) dx
          
          DO IEL=1,IELmax-IELset+1
          
            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
            
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.

              OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
              
              ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
              
              derror = derror + &
                       OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.
          !
          ! It's time to call our coefficient function to calculate the
          ! X-derivative values in the cubature points:  u(x,y)
          ! The result is saved in Dcoefficients(:,:,1)
          
          CALL ffunctionReference (DER_DERIV_X,rdiscretisation, &
                      INT(IELmax-IELset+1),ncubp,p_DcubPtsReal, &
                      IdofsTrial,rintSubset,p_rcollection, &
                      Dcoefficients(:,1:IELmax-IELset+1_I32,1))
                      
          ! Calculate the Y-derivative to Dcoefficients(:,:,2)

          CALL ffunctionReference (DER_DERIV_Y,rdiscretisation, &
                      INT(IELmax-IELset+1),ncubp,p_DcubPtsReal, &
                      IdofsTrial,rintSubset,p_rcollection, &
                      Dcoefficients(:,1:IELmax-IELset+1_I32,2))
          
          ! Calculate the X/Y-derivative of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,3..4)
          
          CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
               p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
               p_elementDistribution%itrialElement, IdofsTrial, &
               ncubp, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV_X,&
               Dcoefficients(:,1:IELmax-IELset+1_I32,3))        

          CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
               p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
               p_elementDistribution%itrialElement, IdofsTrial, &
               ncubp, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV_Y,&
               Dcoefficients(:,1:IELmax-IELset+1_I32,4))        
          
          ! Subtracttion of Dcoefficients(:,:,1..2) from Dcoefficients(:,:,3..4) gives
          ! the error "grad(u-u_h)(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
          
          DO IEL=1,IELmax-IELset+1
          
            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
            
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.

              OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
              
              ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
              
              derror = derror + OM * &
                       ((Dcoefficients(icubp,IEL,3)-Dcoefficients(icubp,IEL,1))**2 + &
                        (Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,2))**2)

            END DO ! ICUBP 

          END DO ! IEL
        
        CASE DEFAULT
          PRINT *,'pperr_scalar2d_conf: Unknown error type identifier!'
          STOP
        END SELECT
    
      END DO ! IELset
      
      ! Release memory
      CALL domint_doneIntegration(rintSubset)

      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial)

    END DO ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    derror = SQRT(derror)

  END SUBROUTINE

END MODULE
