!#########################################################################
!# ***********************************************************************
!# <name> pprocerror </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating errors and norms
!# of finite element functions.
!#
!# The following routines can be found in this module:
!#
!# 1.) pperr_scalar
!#     -> Calculate $L_1$-error, $L_2$-error or $H_1$-error to an
!#        analytic reference function or the $L_1$-norm, $L_2$-norm
!#        or $H_1$-norm of a FE function:
!#   $$ int_\Omega u-u_h dx , \qquad int_\Omega \nabla u-\nabla u_h dx $$
!#
!# 2.) pperr_scalarBoundary2d
!#     -> On a 2D boundary segment, calculate $L_1$-error, $L_2$-error 
!#        or $H_1$-error to an analytic reference function or the 
!#        $L_1$-norm, $L_2$-norm or $H_1$-norm of a FE function.
!#   $$ int_\Gamma u-cu_h dx , \qquad int_\Gamma \nabla u-c\nabla u_h dx $$
!#
!# 3.) pperr_scalarL2ErrorEstimate
!#     -> Calculate $L_2$-error to two different scalar vectors of a 
!#        FE function:
!#   $$ int_\Omega u_h-u_ref dx $$
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is 
!#        some reference solution vector which is supposed to be a 
!#        better approximation of the true solution.
!#
!# 4.) pperr_blockL2ErrorEstimate
!#     -> Calculate $L_2$-error to two different block vectors of a 
!#        FE function:
!#   $$ int_\Omega u_h-u_ref dx $$
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is
!#        some reference solution vector which is supposed to be a 
!#        better approximation of the true solution.
!#
!# 5.) pperr_scalarL1ErrorEstimate
!#     -> Calculate $L_1$-error to two different scalar vectors of a 
!#        FE function:
!#   $$ int_\Omega u_h-u_ref dx $$
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is 
!#        some reference solution vector which is supposed to be a 
!#        better approximation of the true solution.
!#
!# 6.) pperr_blockL1ErrorEstimate
!#     -> Calculate $L_1$-error to two different block vectors of a 
!#        FE function:
!#   $$ int_\Omega u_h-u_ref dx $$
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is
!#        some reference solution vector which is supposed to be a 
!#        better approximation of the true solution.
!#
!# 7.) pperr_scalarStandardDeviation
!#     -> Calculate the standard deviation of a scalar vector
!#
!# 8.) pperr_blockStandardDeviation
!#     -> Calculate the standard deviation of a block vector
!# </purpose>
!#########################################################################

MODULE pprocerror

  USE fsystem
  USE storage
  USE boundary
  USE cubature
  USE triangulation
  USE linearalgebra
  USE linearsystemscalar
  USE linearsystemblock
  USE scalarpde
  USE spatialdiscretisation
  USE domainintegration
  USE collection
  USE elementpreprocessing
  USE feevaluation

  IMPLICIT NONE

!<constants>

!<constantblock description = "Identifiers for the type of error to be computed.">

  ! $L_2$-error/norm
  INTEGER, PARAMETER :: PPERR_L2ERROR = 1
  
  ! $H_1$-error/norm
  INTEGER, PARAMETER :: PPERR_H1ERROR = 2

  ! $L_1$-error/norm
  INTEGER, PARAMETER :: PPERR_L1ERROR = 3
  
  
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
  ! This routine calculates the error or the norm, respectively, of a given 
  ! finite element function in rvector to a given analytical 
  ! callback function ffunctionReference.
  !
  ! If ffunctionReference is specified, the routine calculates
  !   $$ ||y-z||_{L_2}  \textrm{ or }  ||y-z||_{H_1}$$
  ! with $y$=rvectorScalar and $z$=ffunctionReference.
  !
  ! If ffunctionReference is not specified, the routine calculates
  !   $$ ||y||_{L_2}  \textrm{ or }  ||y||_{H_1}.$$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
  !
  ! If the H1-error is desired and the element does not provide first
  ! derivatives, then this routine sets derror to -1 to indicate that the
  ! calculation of the H1-error is not available.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN), TARGET :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_refFunctionSc.inc'
  OPTIONAL :: ffunctionReference
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! If not specified, the discretisation structure in the vector is used.
  ! If specified, the discretisation structure must be 'compatible' to the
  ! vector (concerning NEQ,...). pperr_scalar uses the cubature formula
  ! specifier of the linear form in rdiscretisation to compute the integrals
  ! for the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET, OPTIONAL :: rdiscretisation
!</input>

!<output>
  ! The calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer :: i, imaxder
    integer(I32) :: celement
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! Get the correct discretisation structure and check if we can use it.
    IF (PRESENT(rdiscretisation)) THEN
      p_rdiscretisation => rdiscretisation
      CALL lsyssc_checkDiscretisation (rvectorScalar,p_rdiscretisation)
    ELSE
      p_rdiscretisation => rvectorScalar%p_rspatialdiscr
    END IF
    
    IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
      CALL output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
      CALL sys_halt()
    END IF
    
    ! The vector must be unsorted, otherwise we can't set up the vector.
    IF (rvectorScalar%isortStrategy .GT. 0) THEN
      CALL output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
      CALL sys_halt()
    END IF
    
    ! There is one thing we need to assure: If the user wants to calculate
    ! an H1-error, the corresponding element must have first derivatives!
    if(cerrortype .eq. PPERR_H1ERROR) then
      
      ! So, let's loop through all element distributions of the spatial
      ! discretisation structure.
      do i = 1, p_rdiscretisation%inumFESpaces
        
        ! Get the element of the FE space
        celement = p_rdiscretisation%RelementDistr(i)%celement
        
        ! Get the maximum derivative of the element
        imaxder = elem_getMaxDerivative(celement)
        
        ! If the maximum derivate is 1 then the element's derivatives can
        ! not be evaluated (since they are zero). We will simply return
        ! an error of -1 to indicate that the H1-error cannot be computed.
        if(imaxder .le. 1) then
          derror = -1.0_DP
          return
        end if
      end do
      
      ! If we come out here, then the element provides first derivatives, and
      ! we are ready to calculate the H1-error.
    
    end if
  
    ! Do we have a uniform triangulation? Would simplify a lot...
    IF (p_rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN 
    
      IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
      
        ! Do we have a 1D, 2D or 3D discretisation here?
        SELECT CASE(p_rdiscretisation%ndimension)
        CASE (NDIM1D)
          CALL pperr_scalar1d_conf (rvectorScalar,cerrortype,derror,&
                           p_rdiscretisation,ffunctionReference,rcollection)
        CASE (NDIM2D)
          CALL pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,&
                           p_rdiscretisation,ffunctionReference,rcollection)
        CASE (NDIM3D)
          CALL pperr_scalar3d_conf (rvectorScalar,cerrortype,derror,&
                           p_rdiscretisation,ffunctionReference,rcollection)
        END SELECT
      ELSE
        CALL output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
        CALL sys_halt()
      END IF
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    ELSE IF (p_rdiscretisation%ccomplexity .EQ. SPDISC_CONFORMAL) THEN 
    
      IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
        CALL pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,&
                           p_rdiscretisation,ffunctionReference,rcollection)
      ELSE
        CALL output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
        CALL sys_halt() 
      END IF
    
    ELSE
      CALL output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
      CALL sys_halt()
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalar1d_conf (rvectorScalar,cerrortype,derror,&
                                  rdiscretisation,ffunctionReference,rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 1D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN), TARGET :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_refFunctionSc.inc'
  OPTIONAL :: ffunctionReference
!</input>

!<output>
  ! Array receiving the calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: i,k,icurrentElementDistr, ICUBP, NVE
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
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
    
    ! Current element distribution
    TYPE(t_elementDistribution), POINTER :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    TYPE(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! Type of transformation from the reference to the real element 
    INTEGER :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    INTEGER(I32) :: cevaluationTag

    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial
  
    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .FALSE.
    SELECT CASE (cerrortype)
    CASE (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC1D) = .TRUE.
    CASE (PPERR_H1ERROR) 
      Bder(DER_DERIV1D_X) = .TRUE.
    CASE DEFAULT
      CALL output_line('Unknown error type identifier!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar1d_conf')
      CALL sys_halt()
    END SELECT
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = MIN(PPERR_NELEMSIM,p_rtriangulation%NEL)
                               
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      IF (p_relementDistribution%NEL .EQ. 0) CYCLE

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      DO i=1,ncubp
        DO k=1,UBOUND(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        END DO
      END DO

      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,2))
    
      ! Initialisation of the element set.
      CALL elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      IF (PRESENT(ffunctionReference)) THEN
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REALPOINTS)
      END IF

      ! Make sure that we have determinants.
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        CALL elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        SELECT CASE (cerrortype)
        
        CASE (PPERR_L2ERROR)
          
          ! L2-error uses only the values of the function.
          
          IF (PRESENT(ffunctionReference)) THEN
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_FUNC1D,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp, &
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset, &
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          END IF

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC1D,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))
                  
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
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
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
              
              derror = derror + &
                       OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          IF (PRESENT(ffunctionReference)) THEN
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_FUNC1D,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          END IF

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC1D,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))

          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx
          
          DO IEL=1,IELmax-IELset+1
          
            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
            
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L1-error is:   int_... abs(u-u_h) dx
              
              derror = derror + &
                       OM * ABS(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.

          IF (PRESENT(ffunctionReference)) THEN          
            ! It's time to call our coefficient function to calculate the
            ! X-derivative values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_DERIV1D_X,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
                        
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1:2) = 0.0_DP
          END IF
          
          ! Calculate the X/Y-derivative of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,3)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV1D_X,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))

          ! Subtraction of Dcoefficients(:,:,1..2) from Dcoefficients(:,:,3..4) gives
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
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
              
              derror = derror + OM * &
                       (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

            END DO ! ICUBP 

          END DO ! IEL
        
        CASE DEFAULT
          CALL output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar1d_conf')
          CALL sys_halt()
        END SELECT
        
        ! Release the temporary domain integration structure again
        call domint_doneIntegration (rintSubset)
    
      END DO ! IELset
      
      ! Release memory
      CALL elprep_releaseElementSet(revalElementSet)

      DEALLOCATE(p_DcubPtsRef)
      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial)

    END DO ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    IF ((cerrortype .EQ. PPERR_L2ERROR) .OR.&
        (cerrortype .EQ. PPERR_H1ERROR)) THEN
      derror = SQRT(derror)
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,&
                                  rdiscretisation,ffunctionReference,rcollection)

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
  
  ! A discretisation structure specifying how to compute the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information for callback routines.
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_refFunctionSc.inc'
  OPTIONAL :: ffunctionReference
!</input>

!<output>
  ! Array receiving the calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: i,k,icurrentElementDistr, ICUBP, NVE
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
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
    
    ! Current element distribution
    TYPE(t_elementDistribution), POINTER :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    TYPE(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    INTEGER :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    INTEGER(I32) :: cevaluationTag

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .FALSE.
    SELECT CASE (cerrortype)
    CASE (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC) = .TRUE.
    CASE (PPERR_H1ERROR) 
      Bder(DER_DERIV_X) = .TRUE.
      Bder(DER_DERIV_Y) = .TRUE.
    CASE DEFAULT
      CALL output_line('Unknown error type identifier!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
      CALL sys_halt()
    END SELECT
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = MIN(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      IF (p_relementDistribution%NEL .EQ. 0) CYCLE

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      DO i=1,ncubp
        DO k=1,UBOUND(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        END DO
      END DO
      
      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,4))
    
      ! Initialisation of the element set.
      CALL elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      IF (PRESENT(ffunctionReference)) THEN
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REALPOINTS)
      END IF
                      
      ! Make sure that we have determinants.
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        CALL elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        SELECT CASE (cerrortype)
        
        CASE (PPERR_L2ERROR)
          
          ! L2-error uses only the values of the function.
          
          IF (PRESENT(ffunctionReference)) THEN
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_FUNC,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          END IF

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))

          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
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
            !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
              
              derror = derror + &
                       OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          IF (PRESENT(ffunctionReference)) THEN
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_FUNC,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          END IF

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx
          
          DO IEL=1,IELmax-IELset+1
          
            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
            
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.

              OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
              
              ! L1-error is:   int_... abs(u-u_h) dx
              
              derror = derror + &
                       OM * ABS(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.

          IF (PRESENT(ffunctionReference)) THEN          
            ! It's time to call our coefficient function to calculate the
            ! X-derivative values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_DERIV_X,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
                        
            ! Calculate the Y-derivative to Dcoefficients(:,:,2)

            CALL ffunctionReference (DER_DERIV_Y,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,2),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1:2) = 0.0_DP
          END IF
          
          ! Calculate the X/Y-derivative of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,3..4)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV_X,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,3))

          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV_Y,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,4))

          ! Subtraction of Dcoefficients(:,:,1..2) from Dcoefficients(:,:,3..4) gives
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
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
              
              derror = derror + OM * &
                       ((Dcoefficients(icubp,IEL,3)-Dcoefficients(icubp,IEL,1))**2 + &
                        (Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,2))**2)

            END DO ! ICUBP 

          END DO ! IEL
        
        CASE DEFAULT
          CALL output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
          CALL sys_halt()
        END SELECT
        
        ! Release the temporary domain integration structure again
        call domint_doneIntegration (rintSubset)
    
      END DO ! IELset
      
      ! Release memory
      CALL elprep_releaseElementSet(revalElementSet)

      DEALLOCATE(p_DcubPtsRef)
      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial)

    END DO ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    IF ((cerrortype .EQ. PPERR_L2ERROR) .OR.&
        (cerrortype .EQ. PPERR_H1ERROR)) THEN
      derror = SQRT(derror)
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalar3d_conf (rvectorScalar,cerrortype,derror,&
                                  rdiscretisation,ffunctionReference,rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 3D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN), TARGET :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information for callback routines.
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_refFunctionSc.inc'
  OPTIONAL :: ffunctionReference
!</input>

!<output>
  ! Array receiving the calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: i,k,icurrentElementDistr, ICUBP, NVE
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
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants 
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
    
    ! Current element distribution
    TYPE(t_elementDistribution), POINTER :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    TYPE(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    INTEGER :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    INTEGER(I32) :: cevaluationTag

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .FALSE.
    SELECT CASE (cerrortype)
    CASE (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC3D) = .TRUE.
    CASE (PPERR_H1ERROR) 
      Bder(DER_DERIV3D_X) = .TRUE.
      Bder(DER_DERIV3D_Y) = .TRUE.
      Bder(DER_DERIV3D_Z) = .TRUE.
    CASE DEFAULT
      CALL output_line('Unknown error type identifier!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar3d_conf')
      CALL sys_halt()
    END SELECT
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = MIN(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      IF (p_relementDistribution%NEL .EQ. 0) CYCLE

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      DO i=1,ncubp
        DO k=1,UBOUND(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        END DO
      END DO
      
      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,6))
    
      ! Initialisation of the element set.
      CALL elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      IF (PRESENT(ffunctionReference)) THEN
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REALPOINTS)
      END IF
                      
      ! Make sure that we have determinants.
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        CALL elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        SELECT CASE (cerrortype)
        
        CASE (PPERR_L2ERROR)
          
          ! L2-error uses only the values of the function.
          
          IF (PRESENT(ffunctionReference)) THEN
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y,z)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_FUNC3D,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          END IF

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y,z).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC3D,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
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
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
              
              derror = derror + &
                       OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          IF (PRESENT(ffunctionReference)) THEN
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y,z)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_FUNC3D,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          END IF

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y,z).
          ! Save the result to Dcoefficients(:,:,2)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC3D,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx
          
          DO IEL=1,IELmax-IELset+1
          
            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
            
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L1-error is:   int_... abs(u-u_h) dx
              
              derror = derror + &
                       OM * ABS(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))

            END DO ! ICUBP 

          END DO ! IEL

        CASE (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.

          IF (PRESENT(ffunctionReference)) THEN          
            ! It's time to call our coefficient function to calculate the
            ! X-derivative values in the cubature points:  u(x,y,z)
            ! The result is saved in Dcoefficients(:,:,1)
            CALL ffunctionReference (DER_DERIV3D_X,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
                        
            ! Calculate the Y-derivative to Dcoefficients(:,:,2)
            CALL ffunctionReference (DER_DERIV3D_Y,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,2),rcollection)

            ! Calculate the Z-derivative to Dcoefficients(:,:,3)
            CALL ffunctionReference (DER_DERIV3D_Z,rdiscretisation, &
                        INT(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial,rintSubset, &
                        Dcoefficients(:,1:IELmax-IELset+1_I32,3),rcollection)
          ELSE
            Dcoefficients(:,1:IELmax-IELset+1_I32,1:3) = 0.0_DP
          END IF
          
          ! Calculate the X/Y/Z-derivative of the FE function in the
          ! cubature points: u_h(x,y,z).
          ! Save the result to Dcoefficients(:,:,4..6)
          
          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_X,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,4))

          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_Y,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,5))

          CALL fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_Z,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,6))

          ! Subtraction of Dcoefficients(:,:,1..3) from Dcoefficients(:,:,4..6) gives
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
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!

              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
              
              derror = derror + OM * &
                       ((Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,1))**2 + &
                        (Dcoefficients(icubp,IEL,5)-Dcoefficients(icubp,IEL,2))**2 + &
                        (Dcoefficients(icubp,IEL,6)-Dcoefficients(icubp,IEL,3))**2)

            END DO ! ICUBP 

          END DO ! IEL
        
        CASE DEFAULT
          CALL output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar3d_conf')
          CALL sys_halt()
        END SELECT
    
        ! Release again the temporary domain integration subset
        call domint_doneIntegration (rintSubset)
    
      END DO ! IELset
      
      ! Release memory
      CALL elprep_releaseElementSet(revalElementSet)

      DEALLOCATE(p_DcubPtsRef)
      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial)

    END DO ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    IF ((cerrortype .EQ. PPERR_L2ERROR) .OR.&
        (cerrortype .EQ. PPERR_H1ERROR)) THEN
      derror = SQRT(derror)
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalarBoundary2D (cerrortype,ccubType,&
      derror,rboundaryRegion,rvectorScalar,ffunctionReference,&
      rcollection,rdiscretisation)

!<description>
  ! This routine calculates the error or the norm, respectively, of a given 
  ! finite element function in rvector to a given analytical 
  ! callback function ffunctionReference.
  !
  ! If ffunctionReference is specified, the routine calculates
  !   $$ ||y-z||_{L_2}  \textrm{ or }  ||y-z||_{L_1}  \textrm{ or }  ||y-z||_{H_1}$$
  ! with $y$=rvectorScalar and $z$=ffunctionReference.
  !
  ! If ffunctionReference is not specified, the routine calculates
  !   $$ ||y||_{L_2}  \textrm{ or }  ||y||_{L_1}  \textrm{ or }  ||y||_{H_1}.$$
  !
  ! If the vector rvectorScalar is not specified, it's assumed to be =0.
  !
  ! rboundaryRegion is a t_boundaryRegion object that allows to
  ! specify the boundary region where the error should be computed.
  ! If not specified, the error is computed over the whole boundary.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  INTEGER, INTENT(IN)                      :: ccubType
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  TYPE(t_boundaryRegion), INTENT(IN), OPTIONAL :: rboundaryRegion
  
  ! OPTIONAL: The FE solution vector. Represents a scalar FE function.
  ! If omitted, the function is assumed to be constantly =1.
  TYPE(t_vectorScalar), INTENT(IN), OPTIONAL, TARGET :: rvectorScalar
  
  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  INCLUDE 'intf_refFunctionScBoundary.inc'
  OPTIONAL :: ffunctionReference
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection

  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! Must be specified if rvectorScalar is not specified as this
  ! describes the domain/triangulation/...
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET, OPTIONAL :: rdiscretisation
!</input>

!<output>
  ! The calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    TYPE(t_boundaryRegion) :: rboundaryReg
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    REAL(DP) :: dlocalError
    INTEGER :: ibdc
    
    ! Get the correct discretisation structure and check if we can use it.
    IF (PRESENT(rdiscretisation)) THEN
      p_rdiscretisation => rdiscretisation
      CALL lsyssc_checkDiscretisation (rvectorScalar,p_rdiscretisation)
    ELSE
      p_rdiscretisation => rvectorScalar%p_rspatialdiscr
    END IF
    
    IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
      CALL output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2D')
      CALL sys_halt()
    END IF
    
    IF (p_rdiscretisation%ndimension .NE. NDIM2D) THEN
      CALL output_line('Only 2D discretisations allowed.',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2D')
      CALL sys_halt()
    END IF

    ! The vector must be unsorted, otherwise we can't set up the vector.
    IF (rvectorScalar%isortStrategy .GT. 0) THEN
      CALL output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2D')
      CALL sys_halt()
    END IF
  
    ! If the boundary region is specified, call pperr_scalarBoundary2d_conf
    ! for that boundary region. Otherwise, call pperr_scalarBoundary2d_conf
    ! for all possible boundary regions and sum up the errors.
    IF (PRESENT(rboundaryRegion)) THEN
      CALL pperr_scalarBoundary2d_conf (cerrortype,ccubType,derror,&
        rboundaryRegion,rvectorScalar,ffunctionReference,&
        p_rdiscretisation,rcollection)
    ELSE
      derror = 0.0_DP
      ! Create a boundary region for each boundary component and call
      ! the calculation routine for that.
      DO ibdc=1,boundary_igetNBoundComp(p_rdiscretisation%p_rboundary)
        CALL boundary_createRegion (p_rdiscretisation%p_rboundary, &
            ibdc, 0, rboundaryReg)
        CALL pperr_scalarBoundary2d_conf (cerrortype,ccubType,dlocalError,&
          rboundaryReg,rvectorScalar,ffunctionReference,&
          p_rdiscretisation,rcollection)
        derror = derror + dlocalError
      END DO
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalarBoundary2d_conf (cerrortype,ccubType,derror,&
      rboundaryRegion,rvectorScalar,ffunctionReference,&
      rdiscretisation,rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  INTEGER, INTENT(IN)                      :: cerrortype
  
  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  INTEGER, INTENT(IN)                      :: ccubType

  ! A t_boundaryRegion specifying the boundary region where
  ! to calculate. 
  TYPE(t_boundaryRegion), INTENT(IN), OPTIONAL :: rboundaryRegion

  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN), OPTIONAL, TARGET :: rvectorScalar
  
  ! A discretisation structure specifying how to compute the error.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information for callback routines.
  TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  ! OPTIONAL: A callback function that provides a coefficient in front
  ! of the FE function. If not specified, a value of 1 is assumed.
  INCLUDE 'intf_refFunctionScBoundary.inc'
  OPTIONAL :: ffunctionReference
!</input>

!<output>
  ! Array receiving the calculated error.
  REAL(DP), INTENT(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: IelementOrientation
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), ALLOCATABLE :: Ielements
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DedgePosition
    
    INTEGER :: ibdc,ibdcoffset,iedge,ilocaledge
    INTEGER(PREC_ELEMENTIDX) :: NEL,NELbdc,iel
    INTEGER(I32) :: ctrafoType
    
    ! The triangulation structure - to shorten some things...
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_EDGEIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundary
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtBoundary
    REAL(DP), DIMENSION(:), POINTER :: p_DedgeParameterValue
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoordinates
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement

    ! Arrays for cubature points
    REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dxi2D,Dpoints,DpointsRef
    REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega1D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dvalues
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE) :: Dcoord
    INTEGER :: ncubp,ipoint,ieltype
    INTEGER(I32) :: icoordSystem
    REAL(DP) :: dlen,dpar1,dpar2
    
    ! Arrays for element distributions for every element
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementDistr

    ! List of element distributions in the discretisation structure
    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: p_RelementDistribution

    ! Get some pointers and arrays for quicker access
    p_rtriangulation => rdiscretisation%p_rtriangulation
    p_RelementDistribution => rdiscretisation%RelementDistr
    
    CALL storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    CALL storage_getbase_int (p_rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    CALL storage_getbase_int (p_rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    CALL storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoordinates)
    CALL storage_getbase_double (p_rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    CALL storage_getbase_double (p_rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
        
    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx
    
    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)
    
    ! Position of the boundary component?
    ibdcoffset = p_IboundaryCpIdx(ibdc)
        
    ! In a first step, we figure out the elements on the boundary and their
    ! orientation. Allocate arrays that are large enough to hold
    ! even all elements on the boundary if necessary.
    ALLOCATE(Ielements(NELbdc), IelementOrientation(NELbdc))
    
    ! Allocate an array saving the start- and end-parameter values
    ! of the edges on the boundary.
    ALLOCATE(DedgePosition(2,NELbdc))
    
    ! Loop through the edges on the boundary component ibdc.
    ! If the edge is inside, remember the element number and figure out
    ! the orientation of the edge.
    ! NEL counts the total number of elements in the region.
    NEL = 0
    DO iedge = 1,NELbdc
      IF (boundary_isInRegion(rboundaryRegion,ibdc,&
          p_DedgeParameterValue(iedge))) THEN
        NEL = NEL + 1
        
        ! Element number
        Ielements(NEL) = p_IelementsAtBoundary(iedge)
        
        ! Element orientation; i.e. the local number of the boundary edge 
        DO ilocaledge = 1,UBOUND(p_IedgesAtElement,1)
          IF (p_IedgesAtElement(ilocaledge,p_IelementsAtBoundary(iedge)) .EQ. &
              p_IedgesAtBoundary(iedge)) EXIT
        END DO
        IelementOrientation(NEL) = ilocaledge
        
        ! Save the start parameter value of the edge -- in length
        ! parametrisation.
        dpar1 = p_DvertexParameterValue(iedge)
        
        ! Save the end parameter value. Be careful: The last edge
        ! must be treated differently!
        IF (iedge .NE. NELbdc) THEN
          dpar2 = p_DvertexParameterValue(iedge+1)
        ELSE
          dpar2 = boundary_dgetMaxParVal(&
            rdiscretisation%p_rboundary,ibdc)
        END IF
        
        DedgePosition(1,NEL) = &
          boundary_convertParameter(rdiscretisation%p_rboundary, &
            ibdc, dpar1, rboundaryRegion%cparType, BDR_PAR_LENGTH)
            
        DedgePosition(2,NEL) = &
          boundary_convertParameter(rdiscretisation%p_rboundary, &
            ibdc, dpar2, rboundaryRegion%cparType, BDR_PAR_LENGTH)
         
      END IF
    END DO
    
    ! Get the parameter values of the 1D cubature formula
    ! as well as the number of cubature points ncubp
    CALL cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega1D)
    
    ! Map the 1D cubature points to the edges in 2D.
    ALLOCATE(Dxi2D(ncubp,NDIM2D+1,NEL))
    IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
      ! All elements have the same type. Get the type of the
      ! corresponding coordinate system and transform the coordinates.
      icoordSystem = elem_igetCoordSystem(&
        rdiscretisation%RelementDistr(1)%celement)
      DO iel = 1,NEL
        CALL trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      END DO
    ELSE
      ! The type of the coordinate system may change with every element.
      ! So we may have to switch... celements in the discretisation
      ! structure informs us about the element type.
      CALL storage_getbase_int (rdiscretisation%h_IelementDistr,&
          p_IelementDistr)
      DO iel = 1,NEL
        ieltype = p_RelementDistribution(p_IelementDistr(Ielements(iel)))%celement
        icoordSystem = elem_igetCoordSystem(ieltype)
        CALL trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      END DO
    END IF
    
    ! Transpose the coordinate array such that we get coordinates
    ! we can work with.
    ALLOCATE(DpointsRef(NDIM2D+1,ncubp,NEL))
    DO iel=1,NEL
      DpointsRef(:,:,iel) = TRANSPOSE(Dxi2D(:,:,iel))
    END DO
    
    ! Dxi2D is not needed anymore.
    DEALLOCATE(Dxi2D)
    
    ! If the reference function exists, calculate the coordinates of the
    ! points on world coordinates
    IF (PRESENT(ffunctionReference)) THEN
      
      ! We need the real coordinates of the points.
      ALLOCATE(Dpoints(ncubp,NDIM2D+1,NEL))
      
      IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
        ! All elements with the same transformation
        ctrafoType = elem_igetTrafoType(&
            rdiscretisation%RelementDistr(1)%celement)
        DO iel = 1,NEL

          ! Get the points forming the element
          DO ipoint = 1,UBOUND(p_IverticesAtElement,1)
            Dcoord(1,ipoint) = &
                p_DvertexCoordinates(1,p_IverticesAtElement(ipoint,iel))
            Dcoord(2,ipoint) = &
                p_DvertexCoordinates(2,p_IverticesAtElement(ipoint,iel))
          END DO

          ! Transform the cubature points
          DO ipoint = 1,ncubp
            CALL trafo_calcRealCoords (ctrafoType,Dcoord,&
                DpointsRef(:,ipoint,iel),Dpoints(:,ipoint,iel))
          END DO
        END DO
      ELSE
        ! Transformation can be different for all elements
        DO iel = 1,NEL

          ! Get the points forming the element
          DO ipoint = 1,UBOUND(p_IverticesAtElement,1)
            Dcoord(1,ipoint) = &
                p_DvertexCoordinates(1,p_IverticesAtElement(ipoint,iel))
            Dcoord(2,ipoint) = &
                p_DvertexCoordinates(2,p_IverticesAtElement(ipoint,iel))
          END DO

          ! Transform the cubature points
          DO ipoint = 1,ncubp
            ieltype = p_RelementDistribution(&
                p_IelementDistr(Ielements(iel)))%celement
            ctrafoType = elem_igetTrafoType(ieltype)
            CALL trafo_calcRealCoords (ctrafoType,Dcoord,&
                DpointsRef(:,ipoint,iel),Dpoints(:,ipoint,iel))
          END DO
        END DO
      END IF
    
    END IF    
    
    ! So Dxi2 defines the coordinates on the reference element for all
    ! elements. Generally said, we have to evaluate the elements in these
    ! points now. That can be done by using fevl_evaluate_mult.
    !
    ! Which type of integral is to calculate? H1 or L2 or L1?
    SELECT CASE (cerrortype)
    CASE (PPERR_L2ERROR)
    
      ALLOCATE (Dvalues(ncubp,NEL,2))
      Dvalues = 0.0_DP
      
      ! If the FE function exists, evaluate it.
      IF (PRESENT(rvectorScalar)) THEN
        DO iel = 1,NEL
          ! Evaluate on the element, write results to Dvalues
          CALL fevl_evaluate_mult (DER_FUNC, Dvalues(:,iel,1), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
        END DO
      END IF

      ! If the reference function exists, evaluate it.
      IF (PRESENT(ffunctionReference)) THEN
        
        ! Evaluate the reference function on the boundary
        CALL ffunctionReference (DER_FUNC,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,2),rcollection)
            
      END IF
      
      ! Linear combination to get the actual values in the cubature points.
      DO iel = 1,NEL
        DO ipoint = 1,ncubp
          Dvalues(ipoint,iel,1) = Dvalues(ipoint,iel,2)-Dvalues(ipoint,iel,1)
        END DO
      END DO
    
      ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
      ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
      ! element. We are finally able to calculate the integral!
      ! That means, run over all the edges and sum up...
      ! (ok, if rvectorScalar is not specified, we have
      !  -u_h(x,y) in Dvalues1(:,:,1), but as we take the square,
      !  it doesn't matter if we have u_h or -u_h there!)
      
      derror = 0.0_DP
      DO iel = 1,NEL
      
        ! Get the length of the edge. Let's use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit inverval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
      
        DO ipoint = 1,ncubp
          derror = derror + dlen * Domega1D(ipoint) * (Dvalues(ipoint,iel,1)**2)
        END DO
      END DO
      
      DEALLOCATE(Dvalues)

    CASE (PPERR_L1ERROR)
    
      ALLOCATE (Dvalues(ncubp,NEL,2))
      Dvalues = 0.0_DP
      
      ! If the FE function exists, evaluate it.
      IF (PRESENT(rvectorScalar)) THEN
        DO iel = 1,NEL
          ! Evaluate on the element, write results to Dvalues
          CALL fevl_evaluate_mult (DER_FUNC, Dvalues(:,iel,1), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
        END DO
      END IF

      ! If the reference function exists, evaluate it.
      IF (PRESENT(ffunctionReference)) THEN
        
        ! Evaluate the reference function on the boundary
        CALL ffunctionReference (DER_FUNC,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,2),rcollection)
            
      END IF
      
      ! Linear combination to get the actual values in the cubature points.
      DO iel = 1,NEL
        DO ipoint = 1,ncubp
          Dvalues(ipoint,iel,1) = Dvalues(ipoint,iel,2)-Dvalues(ipoint,iel,1)
        END DO
      END DO
    
      ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
      ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
      ! element. We are finally able to calculate the integral!
      ! That means, run over all the edges and sum up...
      ! (ok, if rvectorScalar is not specified, we have
      !  -u_h(x,y) in Dvalues1(:,:,1), but as we take the square,
      !  it doesn't matter if we have u_h or -u_h there!)
      
      derror = 0.0_DP
      DO iel = 1,NEL
      
        ! Get the length of the edge. Let's use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit inverval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
      
        DO ipoint = 1,ncubp
          derror = derror + dlen * Domega1D(ipoint) * ABS(Dvalues(ipoint,iel,1))
        END DO
      END DO
      
      DEALLOCATE(Dvalues)
      
    CASE (PPERR_H1ERROR)
    
      ALLOCATE (Dvalues(ncubp,NEL,4))
      Dvalues = 0.0_DP
      
      ! If the FE function exists, evaluate it.
      IF (PRESENT(rvectorScalar)) THEN
        DO iel = 1,NEL
          ! Evaluate on the element, write results to Dvalues.
          !
          ! X-derivative
          CALL fevl_evaluate_mult (DER_DERIV_X, Dvalues(:,iel,1), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
              
          ! Y-derivative
          CALL fevl_evaluate_mult (DER_DERIV_Y, Dvalues(:,iel,2), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
        END DO
      END IF

      ! If the reference function exists, evaluate it.
      IF (PRESENT(ffunctionReference)) THEN

        ! Evaluate the reference function on the boundary
        !
        ! X-derivative
        CALL ffunctionReference (DER_DERIV_X,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,3),rcollection)

        ! Y-derivative
        CALL ffunctionReference (DER_DERIV_Y,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,4),rcollection)
            
      END IF

      ! Linear combination to get the actual values in the cubature points.
      ! ||u-u_h||_H1 = int ( grad(u-u_h) * grad(u-u_h) )
      !              ~ sum grad_x(u-u_h)**2 + grad_y(u-u_h)
      DO iel = 1,NEL
        DO ipoint = 1,ncubp
          Dvalues(ipoint,iel,1) = (Dvalues(ipoint,iel,1)-Dvalues(ipoint,iel,3))**2 + &
                                  (Dvalues(ipoint,iel,2)-Dvalues(ipoint,iel,4))**2 
        END DO
      END DO
      
      ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
      ! "grad(u(x,y))-grad(u_h(x,y))" -- in every cubature point on every
      ! element. We are finally able to calculate the integral!
      ! That means, run over all the edges and sum up...
      
      derror = 0.0_DP
      DO iel = 1,NEL
      
        ! Get the length of the edge. Let's use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit inverval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
      
        DO ipoint = 1,ncubp
          derror = derror + dlen * Domega1D(ipoint) * (Dvalues(ipoint,iel,1)**2)
        END DO
      END DO
      
      DEALLOCATE(Dvalues)  
      
    END SELECT

    ! Release memory

    IF (PRESENT(ffunctionReference)) DEALLOCATE(Dpoints)
      
    DEALLOCATE(DedgePosition)
    DEALLOCATE(Ielements, IelementOrientation)
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalarErrorEstimate (rvector,rvectorRef,ctype,derror,&
                                        rdiscretisationRef,relementError)

!<description>
  ! This routine calculates the L2-error of a given FE function in rvector
  ! and a reference vector given in rvectorRef. Both vectors must have the
  ! same number of blocks. As an example, one can think of the consistent
  ! FE gradient and some recovered reference gradient, c.f. ZZ-technique.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution vector
    TYPE(t_vectorScalar), INTENT(IN), TARGET                    :: rvector

    ! FE reference solution vector
    TYPE(t_vectorScalar), INTENT(IN), TARGET                    :: rvectorRef
            
    ! Type of error to compute. A PPERR_xxERROR constant.
    ! PPERR_L2ERROR computes the L2-error, PPERR_L1ERROR the L1-error.
    INTEGER, INTENT(IN) :: ctype
            
    ! OPTIONAL: A discretisation structure specifying how to compute the error.
    ! If not specified, the discretisation structure in the reference gradient 
    ! is used. If specified, the discretisation structure must be 'compatible'
    ! to the two gradient vectors (concerning NEQ,...). pperr_gradient uses the
    ! cubature formula specifier of the linear form in rdiscretisation to 
    ! compute the integrals for the error.
    TYPE(t_spatialDiscretisation), INTENT(IN), TARGET, OPTIONAL :: rdiscretisationRef
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated error on each element.
    TYPE(t_vectorScalar), INTENT(INOUT), OPTIONAL               :: relementError
!</inputoutput>

!<output>
    ! The calculated error.
    REAL(DP), INTENT(OUT)                                       :: derror 
!</output>
!</subroutine>

    ! local variables
    TYPE(t_vectorBlock) :: rvectorBlock,rvectorBlockRef
    TYPE(t_blockDiscretisation) :: rDiscr,rDiscrRef

    ! Create block discretisations with one component
    IF (ASSOCIATED(rvector%p_rspatialdiscr)) THEN
      CALL spdiscr_createBlockDiscrInd(rvector%p_rspatialdiscr, rDiscr)
    ELSE
      CALL output_line('Vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarL2ErrorEstimate')
      CALL sys_halt()
    END IF

    IF (ASSOCIATED(rvectorRef%p_rspatialdiscr)) THEN
      CALL spdiscr_createBlockDiscrInd(rvectorRef%p_rspatialdiscr, rDiscrRef)
    ELSE
      CALL output_line('Reference vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarL2ErrorEstimate')
      CALL sys_halt()
    END IF

    ! Create block vectors with one block
    CALL lsysbl_createVecFromScalar(rvector, rvectorBlock, rDiscr)
    CALL lsysbl_createVecFromScalar(rvectorRef, rvectorBlockRef, rDiscrRef)

    ! Call block version
    CALL pperr_blockErrorEstimate(rvectorBlock, rvectorBlockRef,ctype,&
        derror, rdiscretisationRef,relementError)

    ! Release auxiliary block discretisations
    CALL spdiscr_releaseBlockDiscr(rDiscr)
    CALL spdiscr_releaseBlockDiscr(rDiscrRef)

    ! Release auxiliary block vectors
    CALL lsysbl_releaseVector(rvectorBlock)
    CALL lsysbl_releaseVector(rvectorBlockRef)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_blockErrorEstimate (rvector,rvectorRef,ctype,derror,&
                                       rdiscretisationRef,relementError)

!<description>
  ! This routine calculates the L2-error of a given FE function in rvector
  ! and a reference vector given in rvectorRef. Both vectors must have the
  ! same number of blocks. As an example, one can think of the consistent
  ! FE gradient and some recovered reference gradient, c.f. ZZ-technique.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution block vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET                     :: rvector

    ! FE reference solution block vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET                     :: rvectorRef
            
    ! Type of error to compute. A PPERR_xxERROR constant.
    ! PPERR_L2ERROR computes the L2-error, PPERR_L1ERROR the L1-error.
    INTEGER, INTENT(IN) :: ctype

    ! OPTIONAL: A discretisation structure specifying how to compute the error.
    ! If not specified, the discretisation structure in the reference gradient 
    ! is used. If specified, the discretisation structure must be 'compatible'
    ! to the two gradient vectors (concerning NEQ,...). pperr_gradient uses the
    ! cubature formula specifier of the linear form in rdiscretisation to 
    ! compute the integrals for the error.
    TYPE(t_spatialDiscretisation), INTENT(IN), TARGET, OPTIONAL :: rdiscretisationRef
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated error on each element.
    TYPE(t_vectorScalar), INTENT(INOUT), OPTIONAL               :: relementError
!</inputoutput>

!<output>
    ! The calculated error.
    REAL(DP), INTENT(OUT)                                       :: derror 
!</output>
!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisationRef
    INTEGER      :: i,k,icurrentElementDistr,iblock,ICUBP,NVE
    INTEGER(I32) :: IEL, IELmax, IELset,IELGlobal
    REAL(DP)     :: OM,delementError

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
    INTEGER :: indofTrial,indofTrialRef
    
    ! The triangulation structure - to shorten some things...
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! A pointer to an element-number list
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcubPtsRef
    
    ! Jacobian determinant in all points
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj

    ! Pointer to the element error
    REAL(DP), DIMENSION(:), POINTER :: p_DelementError
    
    ! Current element distribution
    TYPE(t_elementDistribution), POINTER :: p_relementDistribution,p_relementDistributionRef
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
    
    ! Type of transformation from the reference to the real element 
    INTEGER :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    INTEGER(I32) :: cevaluationTag

    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    TYPE(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial,IdofsTrialRef
    
    ! Get the correct discretisation structure for the solution vector
    p_rdiscretisation => rvector%p_rblockDiscr%RspatialDiscr(1)
    DO iblock=2,rvector%nblocks
      CALL lsyssc_checkDiscretisation (rvector%RvectorBlock(iblock),p_rdiscretisation)
    END DO

    ! Get the correct discretisation structure for the reference 
    ! vector and check if we can use it.
    IF (PRESENT(rdiscretisationRef)) THEN
      p_rdiscretisationRef => rdiscretisationRef
      CALL lsyssc_checkDiscretisation (rvectorRef%RvectorBlock(1),p_rdiscretisationRef)
    ELSE
      p_rdiscretisationRef => rvectorRef%p_rblockDiscr%RspatialDiscr(1)
    END IF
    DO iblock=2,rvectorRef%nblocks
      CALL lsyssc_checkDiscretisation (rvectorRef%RvectorBlock(iblock),p_rdiscretisationRef)
    END DO
    
    IF (.NOT. ASSOCIATED(p_rdiscretisation) .OR.&
        .NOT. ASSOCIATED(p_rdiscretisationRef)) THEN
      CALL output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
      CALL sys_halt()
    END IF
    
    ! The vectors must have the same number of blocks
    IF (rvector%nblocks .NE. rvectorRef%nblocks) THEN
      CALL output_line('Vectors have different number of blocks!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
      CALL sys_halt()
    END IF
    
    ! The vector must be unsorted.
    DO iblock=1,rvector%nblocks
      IF (rvector%RvectorBlock(iblock)%isortStrategy    .GT. 0 .OR.&
          rvectorRef%RvectorBlock(iblock)%isortStrategy .GT. 0) THEN
        CALL output_line('Vectors must be unsorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
        CALL sys_halt()
      END IF
    END DO

    ! We only need the function values of basis functions
    Bder = .FALSE.
    Bder(DER_FUNC) = .TRUE.

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscretisationRef%p_rtriangulation

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = MIN(PPERR_NELEMSIM,p_rtriangulation%NEL)

    ! Set the current error to 0 and add the error contributions of each element to that.
    derror = 0.0_DP

    ! Get a pointer to the element error (if required)
    IF (PRESENT(relementError)) THEN
      CALL lsyssc_getbase_double(relementError,p_DelementError)
      CALL lalg_clearVectorDble (p_DelementError)
    END IF

    ! Check that both discretisations have the same number of element distributions
    IF (p_rdiscretisation%inumFESpaces .NE. &
        p_rdiscretisationRef%inumFESpaces) THEN
      CALL output_line('Number of element distributions mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
      CALL sys_halt()
    END IF

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution    => p_rdiscretisation%RelementDistr(icurrentElementDistr)
      p_relementDistributionRef => p_rdiscretisationRef%RelementDistr(icurrentElementDistr)

      ! Check if element distrbutions have different number of elements
      IF (p_relementDistribution%NEL .NE. &
          p_relementDistributionRef%NEL) THEN
        CALL output_line('Number of elements in distributions mismatch!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
        CALL sys_halt()
      END IF

      ! Cancel if this element distribution is empty.
      IF (p_relementDistribution%NEL .EQ. 0) CYCLE

      ! Get the number of local DOF's for trial functions
      indofTrial    = elem_igetNDofLoc(p_relementDistribution%celement)
      indofTrialRef = elem_igetNDofLoc(p_relementDistributionRef%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
     
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_relementDistributionRef%ccubTypeEval, ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistributionRef%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      DO i=1,ncubp
        DO k=1,UBOUND(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        END DO
      END DO
      
      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
      ALLOCATE(IdofsTrialRef(indofTrialRef,nelementsPerBlock))

      ! Allocate memory for the coefficients, that is, two times the spatial dimension
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,2*rvector%nblocks))
    
      ! Initialisation of the element set.
      CALL elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      cevaluationTag = IOR(cevaluationTag,&
                      elem_getEvaluationTag(p_relementDistributionRef%celement))

      ! Make sure that we have determinants.
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)
    
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPERR_NELEMSIM
  
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
        CALL dof_locGlobMapping_mult(p_rdiscretisationRef, p_IelementList(IELset:IELmax), &
                                     IdofsTrialRef)
                                     
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        CALL elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj
            
        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))

        ! Depending on ctype, choose now the error to compute.
        SELECT CASE (ctype)
        CASE (PPERR_L1ERROR)

          ! L1-error uses only the values of the function.

          ! Calculate the values of the FE solution vector and the reference solution 
          ! vector in the cubature points: u_h(x,y) and u_ref(x,y)
          ! Save the result to Dcoefficients(:,:,2*iblock-1) and 
          ! Dcoefficients(:,:,2*iblock)

          DO iblock=1,rvector%nblocks
            
            ! solution vector
            CALL fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1_I32,2*iblock))

            ! solution reference vector
            CALL fevl_evaluate_sim3 (rvectorRef%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistributionRef%celement, IdofsTrialRef, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1_I32,2*iblock-1))

          END DO

          ! Subtraction of Dcoefficients(:,:,2*iblock-1) from Dcoefficients(:,:,2*iblock)
          ! and summing over all iblock=1,..,nblocks gives the error 
          ! $u_h(cubature pt.) - u_ref(cubature pt.)$
          
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u_h-u_ref,u_h-u_ref) dx

          DO IEL=1,IELmax-IELset+1

            ! Initialise element error by 0
            delementError = 0.0_DP

            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
              
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              
              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L1-error is:   int_... abs(u_h-u_ref) dx
              
              DO iblock=1,rvector%nblocks
                delementError = delementError + &
                    OM * ABS(Dcoefficients(icubp,IEL,2*iblock-1)-&
                            Dcoefficients(icubp,IEL,2*iblock))
              END DO

            END DO ! ICUBP 

            ! Apply to global error
            derror = derror + delementError

            ! Store in element error (if required)
            IF (PRESENT(relementError)) THEN
              IELGlobal = p_IelementList(IELset+IEL-1)
              p_DelementError(IELGlobal) = delementError
            END IF

          END DO ! IEL
        
        CASE (PPERR_L2ERROR)
        
          ! L2-error uses only the values of the function.

          ! Calculate the values of the FE solution vector and the reference solution 
          ! vector in the cubature points: u_h(x,y) and u_ref(x,y)
          ! Save the result to Dcoefficients(:,:,2*iblock-1) and 
          ! Dcoefficients(:,:,2*iblock)

          DO iblock=1,rvector%nblocks
            
            ! solution vector
            CALL fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1_I32,2*iblock))

            ! solution reference vector
            CALL fevl_evaluate_sim3 (rvectorRef%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistributionRef%celement, IdofsTrialRef, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1_I32,2*iblock-1))

          END DO

          ! Subtraction of Dcoefficients(:,:,2*iblock-1) from Dcoefficients(:,:,2*iblock)
          ! and summing over all iblock=1,..,nblocks gives the error 
          ! $u_h(cubature pt.) - u_ref(cubature pt.)$
          
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u_h-u_ref,u_h-u_ref) dx

          DO IEL=1,IELmax-IELset+1

            ! Initialise element error by 0
            delementError = 0.0_DP

            ! Loop over all cubature points on the current element
            DO icubp = 1, ncubp
              
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              
              OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
              
              ! L2-error is:   int_... (u_h-u_ref)*(u_h-u_ref) dx
              
              DO iblock=1,rvector%nblocks
                delementError = delementError + &
                    OM * (Dcoefficients(icubp,IEL,2*iblock-1)-&
                          Dcoefficients(icubp,IEL,2*iblock))**2
              END DO

            END DO ! ICUBP 

            ! Apply to global error
            derror = derror + delementError

            ! Store in element error (if required)
            IF (PRESENT(relementError)) THEN
              IELGlobal = p_IelementList(IELset+IEL-1)
              p_DelementError(IELGlobal) = SQRT(delementError)
            END IF

          END DO ! IEL
          
        CASE DEFAULT
          
          CALL output_line('Requested error estimate not implemented!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
          CALL sys_halt()

        END SELECT
        
      END DO ! IELset
      
      ! Release memory
      CALL elprep_releaseElementSet(revalElementSet)

      DEALLOCATE(p_DcubPtsRef)
      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial,IdofsTrialRef)
      
    END DO ! icurrentElementDistr

    IF (ctype .NE. PPERR_L1ERROR) THEN
      ! derror is ||error||^2, so take the square root at last.
      derror = SQRT(derror)
    END IF

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_scalarStandardDeviation (rvector,ddeviation,relementDeviation)

!<description>
  ! This routine calculates the standard deviation
  ! 
  ! $$ \sigma=\sqrt{\int_\Omega r^2 u dx} $$
  !
  ! of a given FE function $u$ in rvector, whereby
  !
  ! $$ r^2=(x-\hat x)^2 + (y-\hat y)^2 + (z-\hat z)^2 $$
  !
  ! and each component is computed from the following relation
  !
  ! $$ \hat x=\int_\Omega x u dx $$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution vector
    TYPE(t_vectorScalar), INTENT(IN), TARGET                    :: rvector
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated deviation on each element.
    TYPE(t_vectorScalar), INTENT(INOUT), OPTIONAL               :: relementDeviation
!</inputoutput>

!<output>
    ! The calculated standard deviation.
    REAL(DP), INTENT(OUT)                                       :: ddeviation
!</output>
!</subroutine>

    ! local variables
    TYPE(t_vectorBlock) :: rvectorBlock
    TYPE(t_blockDiscretisation) :: rDiscr

    ! Create block discretisations with one component
    IF (ASSOCIATED(rvector%p_rspatialdiscr)) THEN
      CALL spdiscr_createBlockDiscrInd(rvector%p_rspatialdiscr, rDiscr)
    ELSE
      CALL output_line('Vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarStandardDeviation')
      CALL sys_halt()
    END IF

    ! Create block vectors with one block
    CALL lsysbl_createVecFromScalar(rvector, rvectorBlock, rDiscr)

    ! Call block version
    CALL pperr_blockStandardDeviation(rvectorBlock, ddeviation, relementDeviation)

    ! Release auxiliary block discretisations
    CALL spdiscr_releaseBlockDiscr(rDiscr)

    ! Release auxiliary block vectors
    CALL lsysbl_releaseVector(rvectorBlock)

  END SUBROUTINE pperr_scalarStandardDeviation

  !****************************************************************************

!<subroutine>

  SUBROUTINE pperr_blockStandardDeviation (rvector,ddeviation,relementDeviation)

!<description>
  ! This routine calculates the standard deviation
  ! 
  ! $$ \sigma=\sqrt{\int_\Omega r^2 u dx} $$
  !
  ! of a given FE function $u$ in rvector, whereby
  !
  ! $$ r^2=(x-\hat x)^2 + (y-\hat y)^2 + (z-\hat z)^2 $$
  !
  ! and each component is computed from the following relation
  !
  ! $$ \hat x=\int_\Omega x u dx $$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution block vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET                     :: rvector
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated deviation on each element.
    TYPE(t_vectorScalar), INTENT(INOUT), OPTIONAL               :: relementDeviation
!</inputoutput>

!<output>
    ! The calculated deviation.
    REAL(DP), INTENT(OUT)                                       :: ddeviation
!</output>
!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    INTEGER      :: i,k,icurrentElementDistr,iblock,ICUBP,NVE,idim
    INTEGER(I32) :: IEL, IELmax, IELset,IELGlobal
    REAL(DP)     :: OM,delementDeviation

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
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcubPtsRef

    ! An array receiving the coordinates of cubature points on
    ! the real element for all elements in a set.
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal

    ! Arrays for saving Jacobian determinants and matrices
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
    
    ! Pointer to the element deviation
    REAL(DP), DIMENSION(:), POINTER :: p_DelementDeviation
    
    ! Current element distribution
    TYPE(t_elementDistribution), POINTER :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients

    ! Mathematical expectation of the center of mass
    REAL(DP), DIMENSION(NDIM3D) :: DmassCenter
    
    ! Type of transformation from the reference to the real element 
    INTEGER :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    INTEGER(I32) :: cevaluationTag
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    TYPE(t_evalElementSet) :: revalElementSet

    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial
    
    ! Get the correct discretisation structure for the solution vector
    p_rdiscretisation => rvector%p_rblockDiscr%RspatialDiscr(1)
    DO iblock=2,rvector%nblocks
      CALL lsyssc_checkDiscretisation (rvector%RvectorBlock(iblock),p_rdiscretisation)
    END DO

    IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
      CALL output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockStandardDeviation')
      CALL sys_halt()
    END IF

    ! The vector must be unsorted.
    DO iblock=1,rvector%nblocks
      IF (rvector%RvectorBlock(iblock)%isortStrategy .GT. 0) THEN
        CALL output_line('Vectors must be unsorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockStandardDeviation')
        CALL sys_halt()
      END IF
    END DO

    ! We only need the function values of basis functions
    Bder = .FALSE.
    Bder(DER_FUNC) = .TRUE.

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscretisation%p_rtriangulation

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = MIN(PPERR_NELEMSIM,p_rtriangulation%NEL)

    ! Set the mathematical expectation of the center of mass to 0
    DmassCenter = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => p_rdiscretisation%RelementDistr(icurrentElementDistr)

      ! Cancel if this element distribution is empty.
      IF (p_relementDistribution%NEL .EQ. 0) CYCLE

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
     
      ! Initialise the cubature formula.
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_relementDistribution%ccubTypeEval,&
          ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      DO i=1,ncubp
        DO k=1,UBOUND(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        END DO
      END DO

      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients.
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,rvector%nblocks))
    
      ! Initialisation of the element set.
      CALL elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! Make sure that we have determinants.
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPERR_NELEMSIM
  
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        CALL elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj
        p_DcubPtsReal => revalElementSet%p_DpointsReal
            
        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))

        ! Standard deviation uses only the values of the function.

        ! Calculate the values of the FE solution vector in the cubature 
        ! points u_h(x,y) and save the result to Dcoefficients(:,:,iblock)

        DO iblock=1,rvector%nblocks
          
          ! solution vector
          CALL fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,iblock))
                  
        END DO

        ! Calculate the mathematical expectation of the center of mass
        ! $\hat x_h=\int_\Omega x u_h dx$

        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega x*u_h dx, for x,y and z
        
        DO IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          DO icubp = 1, ncubp
            
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            
            OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
            
            ! Mathematical expectation of the center of mass is:
            ! int_... x*u_h dx
            
            DO iblock=1,rvector%nblocks
              DO idim=1,p_rdiscretisation%ndimension
                DmassCenter(idim) = DmassCenter(idim) + &
                    OM * p_DcubPtsReal(idim,icubp,IEL) * &
                         Dcoefficients(icubp,IEL,iblock)
              END DO
            END DO

          END DO ! ICUBP 
        
        END DO ! IEL
        
      END DO ! IELset
      
      ! Release memory
      CALL elprep_releaseElementSet(revalElementSet)

      DEALLOCATE(p_DcubPtsRef)      
      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial)
      
    END DO ! icurrentElementDistr

    ! Ok, we have the mathematical expectation of the center of mass.
    ! Let's compute the standard deviation.

    Ddeviation = 0.0_DP

    ! Get a pointer to the element deviation (if required)
    IF (PRESENT(relementDeviation)) THEN
      CALL lsyssc_getbase_double(relementDeviation,p_DelementDeviation)
      CALL lalg_clearVectorDble (p_DelementDeviation)
    END IF

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Make sure that we have determinants.
    cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution    => p_rdiscretisation%RelementDistr(icurrentElementDistr)

      ! Cancel if this element distribution is empty.
      IF (p_relementDistribution%NEL .EQ. 0) CYCLE

      ! Get the number of local DOF's for trial functions
      indofTrial    = elem_igetNDofLoc(p_relementDistribution%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
     
      ! Initialise the cubature formula.
      ! Get cubature weights and point coordinates on the reference element
      CALL cub_getCubPoints(p_relementDistribution%ccubTypeEval,&
          ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      DO i=1,ncubp
        DO k=1,UBOUND(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        END DO
      END DO

      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients.
      ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,rvector%nblocks))
    
      ! Initialisation of the element set.
      CALL elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! Make sure that we have determinants.
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPERR_NELEMSIM
  
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        CALL elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj
        p_DcubPtsReal => revalElementSet%p_DpointsReal

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))

        ! Standard deviation uses only the values of the function.

        ! Calculate the values of the FE solution vector in the cubature 
        ! points u_h(x,y) and save the result to Dcoefficients(:,:,iblock)

        DO iblock=1,rvector%nblocks
          
          ! solution vector
          CALL fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,iblock))
                  
        END DO

        ! Calculate the standard deviation
        ! $\int_\Omega ((x-\hat x_h)^2 + (y-\hat y_h)^2 + (z-\hat z_h)^2) u_h dx$

        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega (x-\hat x_h)*(x-\hat x_h)*u_h dx
        ! and sum up for x,y and z
        
        DO IEL=1,IELmax-IELset+1

          ! Initialise element deviation by 0
          delementDeviation = 0.0_DP

          ! Loop over all cubature points on the current element
          DO icubp = 1, ncubp
            
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            
            OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
            
            ! Standard deviation is: int_... (x-\hat x_h)*(x-\hat x_h)*u_h dx
            ! summed up for all x,y and z
            
            DO iblock=1,rvector%nblocks
              DO idim=1,p_rdiscretisation%ndimension
                delementDeviation = delementDeviation + &
                    OM * Dcoefficients(icubp,IEL,iblock) * &
                        (p_DcubPtsReal(idim,icubp,IEL)-DmassCenter(idim))**2
              END DO
            END DO

          END DO ! ICUBP 
        
          ! Apply to global deviation
          ddeviation = ddeviation + delementDeviation

          ! Store in element deviation (if required)
          IF (PRESENT(relementDeviation)) THEN
            IELGlobal = p_IelementList(IELset+IEL-1)
            p_DelementDeviation(IELGlobal) = SQRT(delementDeviation)
          END IF

        END DO ! IEL
        
      END DO ! IELset
      
      ! Release memory
      CALL elprep_releaseElementSet(revalElementSet)

      DEALLOCATE(p_DcubPtsRef)      
      DEALLOCATE(Dcoefficients)
      DEALLOCATE(IdofsTrial)
      
    END DO ! icurrentElementDistr
    
    ! ddeviation is ||deviation||^2, so take the square root at last.
    IF (ddeviation .GE. HUGE(ddeviation)) THEN
      ddeviation = 0._DP
    ELSE
      ddeviation = SQRT(ddeviation)
    END IF

  END SUBROUTINE pperr_blockStandardDeviation
END MODULE
