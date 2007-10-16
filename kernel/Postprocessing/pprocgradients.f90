!#########################################################################
!# ***********************************************************************
!# <name> pprocgradients </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating gradients
!# of finite element functions.
!#
!# The following routines can be found in this module:
!#
!# 1.) ppgrd_calcGradient
!#     -> Calculates the recovered gradient of a scalar finite element
!#        function.
!#
!# Auxiliary routines, called internally.
!#
!# 1.) ppgrd_calcGrad2DInterpP12Q12cnf
!#     -> Calculate the reconstructed gradient as P1, Q1, P2 or Q2
!#        vector for an arbitrary conformal discretisation, 2D.
!#        Uses 1st order interpolation to reconstruct the gradient.
!# </purpose>
!#########################################################################

MODULE pprocgradients

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
  USE genoutput

  IMPLICIT NONE

!<constants>

!<constantblock description = "Identifiers for the method how to calculate a gradient vector.">

  ! Use standard interpolation to calculate a gradient vector. 1st order. 
  INTEGER, PARAMETER :: PPGRD_INTERPOL = 0
  
  ! ZZ-technique for recovering a gradient. Only usable for special-type
  ! meshes, e.g. pure quad meshes. 2nd order on regular meshes.
  INTEGER, PARAMETER :: PPGRD_ZZTECHNIQUE = 1
  
!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! Number of elements to handle simultaneously when building vectors
  INTEGER :: PPGRD_NELEMSIM   = 1000
  
!</constantblock>

!</constants>

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppgrd_calcGradient (rvectorScalar,rvectorGradient,cgradType)

!<description>
  ! Calculates the recovered gradient of a scalar finite element function.
  ! cgradType decides about the method to use for the calculation.
  ! This parameter is optional; if not specified, a standard method
  ! will be taken.
  !
  ! rvectorGradient receives the reconstructed gradient. For a 2D discretisation,
  ! this must be a 2D vector. For a 3D discretisation, this must be a 3D vector.
  ! The vector must provide a discretisation structure that defines the
  ! finite element space the reconstructed gradient should be calculated in.
  !
  ! Note: Currently, rvectorGradient must be prepared as $P_0,P_1$ or $Q_0,Q_1$
  ! vector, respectively, other types of destination vectors are not allowed! 
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN)         :: rvectorScalar
  
  ! OPTIONAL: Identifier for the method to use for calculating the gradient.
  ! One of the PPGRD_xxxx constants.
  ! If not specified, PPGRD_INTERPOL is taken as default.
  INTEGER, INTENT(IN), OPTIONAL            :: cgradType
!</input>

!<inputoutput>
  ! A block vector receiving the gradient.
  ! The first subvector receives the X-gradient.
  ! In 2D/3D discretisations, the 2nd subvector recevies the Y-gradient.
  ! In 3D discretisations, the 3rd subvector receives the Z-gradient.
  ! The vector must be prepared with a discretisation structure that defines
  ! the destination finite element space for the gradient field.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvectorGradient
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i,j
    LOGICAL :: bisP0,bisQ0,bisQ1, bisP1, bisQ2, bisP2, bisDifferent
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscr
    INTEGER :: imethod

    ! Some basic checks:
    
    imethod = PPGRD_INTERPOL
    IF (PRESENT(cgradType)) imethod = cgradType
    
    ! Dimension of triangulation must be less or equal than number of subvectors
    ! in the block vector.
    
    IF (.NOT. ASSOCIATED(rvectorScalar%p_rspatialDiscretisation)) THEN
      CALL output_line ('No discretisation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      CALL sys_halt()
    END IF
    
    IF ((rvectorScalar%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND.&
        (rvectorScalar%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
      CALL output_line ('Only uniform and conformal discretisations supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      CALL sys_halt()
    END IF
    
    IF (.NOT. ASSOCIATED(rvectorScalar%p_rspatialDiscretisation%p_rtriangulation)) THEN
      CALL output_line ('No triangulation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      CALL sys_halt()
    END IF
    
    IF (rvectorScalar%p_rspatialDiscretisation%p_rtriangulation%ndim .GT. &
        rvectorGradient%nblocks) THEN
      CALL output_line ('Dimension of destination vector not large enough!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      CALL sys_halt()
    END IF
    
    ! There must be given discretisation structures in the destination vector.
    IF (.NOT. ASSOCIATED(rvectorScalar%p_rspatialDiscretisation)) THEN
      CALL output_line ('No discretisation attached to the destination vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      CALL sys_halt()
    END IF

    SELECT CASE (rvectorScalar%p_rspatialDiscretisation%p_rtriangulation%ndim)
    CASE (NDIM2D)
      ! Currently, the destinatino vector must be either pure Q1 or pure P1 or
      ! mixed Q1/P1 -- everything else is currently not supported.
      bisQ0 = .FALSE.
      bisP0 = .FALSE.
      bisQ1 = .FALSE.
      bisP1 = .FALSE.
      bisQ2 = .FALSE.
      bisP2 = .FALSE.
      bisDifferent = .FALSE.
      ! We only check the first two subvectors which we overwrite.
      ! Any additional subvectors are ignored!
      DO i=1,MIN(2,rvectorGradient%nblocks)
        p_rdiscr => rvectorGradient%p_rblockDiscretisation%RspatialDiscretisation(i)
        DO j=1,p_rdiscr%inumFESpaces
          SELECT CASE (&
              elem_getPrimaryElement (p_rdiscr%RelementDistribution(j)%itrialElement))
          CASE (EL_Q0)
            bisQ0 = .TRUE.
          CASE (EL_P0)
            bisP0 = .TRUE.
          CASE (EL_Q1)
            bisQ1 = .TRUE.
          CASE (EL_P1)
            bisP1 = .TRUE.
          CASE (EL_Q2)
            bisQ2 = .TRUE.
          CASE (EL_P2)
            bisP2 = .TRUE.
          CASE DEFAULT
            bisDifferent = .TRUE.
          END SELECT
        END DO
      END DO
      
      IF (bisDifferent) THEN
        CALL output_line ('Only Q0, Q1, P0, P1, Q2 and P2 supported as&
            & discretisation for the destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
        CALL sys_halt()
      END IF
      
      ! The bisXXXX flags might get important later...
      
      ! Depending on the method choosed, call the appropriate gradient
      ! recovery routine.
      SELECT CASE (imethod)
      CASE (PPGRD_INTERPOL)
        CALL ppgrd_calcGrad2DInterpP12Q12cnf (rvectorScalar,rvectorGradient)
      CASE (PPGRD_ZZTECHNIQUE)
        CALL output_line ('ZZ not implemented!!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
        CALL sys_halt()
      END SELECT

    CASE DEFAULT
      CALL output_line ('Unsupported dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      CALL sys_halt()
    END SELECT

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppgrd_calcGrad2DInterpP12Q12cnf (rvector,rvectorGradient)

!<description>
  ! Calculates the recvered gradient of a scalar finite element function
  ! by standard interpolation. Supports conformal 2D discretisations
  ! with $P_1$, $Q_1$, $P_2$ and $Q_2$ mixed in the destination vector.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  TYPE(t_vectorScalar), INTENT(IN)         :: rvector
!</input>

!<inputoutput>
  ! A block vector receiving the gradient.
  ! The first subvector receives the X-gradient.
  ! In 2D/3D discretisations, the 2nd subvector recevies the Y-gradient.
  ! In 3D discretisations, the 3rd subvector receives the Z-gradient.
  ! The vector must be prepared with a discretisation structure that defines
  ! the destination finite element space for the gradient field.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvectorGradient
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i,j,k,icurrentElementDistr, NVE
    LOGICAL :: bnonparTrial
    INTEGER(I32) :: IELmax, IELset
    
    ! Array to tell the element which derivatives to calculate
    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
    
    ! Cubature point coordinates on the reference element
    REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! Cubature formula weights. The cotent is actually not used here.
    REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
    
    ! Number of 'cubature points'. As we acutally don't do cubature
    ! here, this coincides with the number of DOF's on each element
    ! in the destination space.
    INTEGER :: nlocalDOFsDest
    
    ! Number of local degees of freedom for test functions
    INTEGER :: indofTrial,indofDest
    
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
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    
    ! Current element distribution in source- and destination vector
    TYPE(t_elementDistribution), POINTER :: p_elementDistribution
    TYPE(t_elementDistribution), POINTER :: p_elementDistrDest
    
    ! Pointer to the values of the function that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dderivatives
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    INTEGER :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    TYPE(t_domainIntSubset) :: rintSubset
    
    ! An allocateable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrial
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsDest
    
    ! Pointers to the X- and Y-derivative vector
    REAL(DP), DIMENSION(:), POINTER :: p_DxDeriv, p_DyDeriv
    
    ! Pointer to an array that counts the number of elements adjacent to a vertex.
    ! Ok, there's the same information in the triangulation, but that's not
    ! based on DOF's! Actually, we'll calculate how often we touch each DOF 
    ! in the destination space.
    INTEGER :: h_IcontributionsAtDOF
    INTEGER(I32), DIMENSION(:), POINTER :: p_IcontributionsAtDOF
    
    ! Discretisation structures for the source- and destination vector(s)
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrSource, p_rdiscrDest
    
    ! Evaluate the first derivative of the FE functions.

    Bder = .FALSE.
    Bder(DER_DERIV_X) = .TRUE.
    Bder(DER_DERIV_Y) = .TRUE.
    
    ! Get the discretisation structures of the source- and destination space.
    ! Note that we assume here that the X- and Y-derivative is discretised
    ! the same way!
    p_rdiscrSource => rvector%p_rspatialDiscretisation
    p_rdiscrDest => rvectorGradient%p_rblockDiscretisation%RspatialDiscretisation(1)
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscrSource%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = MIN(PPGRD_NELEMSIM,p_rtriangulation%NEL)
    
    ! Get a pointer to the KVERT and DCORVG array
    CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                               p_IverticesAtElement)
    CALL storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, &
                               p_DvertexCoords)
                               
    ! Get pointetrs to the X- and Y-derivative destination vector.
    CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
    CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
    
    ! Array that allows the calculation about the number of elements
    ! meeting in a vertex, based onm DOF's.
    CALL storage_new ('ppgrd_calcGrad2DInterpP1Q1cnf','DOFContrAtVertex',&
                      dof_igetNDofGlob(p_rdiscrDest),&
                      ST_INT, h_IcontributionsAtDOF, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int (h_IcontributionsAtDOF,p_IcontributionsAtDOF)
    
    ! Clear the destination vectors.
    CALL lalg_clearVectorDble (p_DxDeriv)
    CALL lalg_clearVectorDble (p_DyDeriv)
                               
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,p_rdiscrSource%inumFESpaces
    
      ! Activate the current element distribution
      p_elementDistribution => p_rdiscrSource%RelementDistribution(icurrentElementDistr)
      p_elementDistrDest => p_rdiscrDest%RelementDistribution(icurrentElementDistr)
    
      ! Get the number of local DOF's for trial functions
      ! in the source and destination vector.
      indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
      indofDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_elementDistribution%itrialElement)
      
      ! Initialise the cubature formula,
      ! That's a speciat trick here! The FE space of the destination vector
      ! is either P1 or Q1. We create the gradients in the corners of these points
      ! by taking the mean of the gradients of the source vector!
      !
      ! Gor this purpose, we initialise the 'trapezoidal rule' as cubature
      ! formula. Ok, we are not using cubature at all here (thus ignoring
      ! the cubature weights completely), but that way we get the
      ! coordinates of the corners on the reference element automatically --
      ! which coincide with the points where we want to create the gradients!
      !
      ! Note: The returned nlocalDOFsDest will coincide with the number of local DOF's
      ! on each element indofDest!
      SELECT CASE (elem_getPrimaryElement(p_elementDistrDest%itrialElement))
      CASE (EL_P0)
        CALL cub_getCubPoints(CUB_G1_T, nlocalDOFsDest, Dxi, Domega)
      CASE (EL_Q0)
        CALL cub_getCubPoints(CUB_G1X1, nlocalDOFsDest, Dxi, Domega)
      CASE (EL_P1)
        CALL cub_getCubPoints(CUB_TRZ_T, nlocalDOFsDest, Dxi, Domega)
      CASE (EL_Q1)
        CALL cub_getCubPoints(CUB_TRZ, nlocalDOFsDest, Dxi, Domega)
      CASE (EL_P2)
        ! Manually calculate the coordinates of the corners/midpoints on
        ! the reference element.
        Dxi(1,1)  =  1.0_DP
        Dxi(1,2)  =  0.0_DP
        Dxi(1,3)  =  0.0_DP

        Dxi(2,1)  =  0.0_DP
        Dxi(2,2)  =  1.0_DP
        Dxi(2,3)  =  0.0_DP
        
        Dxi(3,1)  =  0.0_DP
        Dxi(3,2)  =  0.0_DP
        Dxi(3,3)  =  1.0_DP
        
        Dxi(4,1)  =  0.5_DP
        Dxi(4,2)  =  0.5_DP
        Dxi(4,3)  =  0.0_DP

        Dxi(5,1)  =  0.0_DP
        Dxi(5,2)  =  0.5_DP
        Dxi(5,3)  =  0.5_DP

        Dxi(6,1)  =  0.5_DP
        Dxi(6,2)  =  0.0_DP
        Dxi(6,3)  =  0.5_DP
        
        nlocalDOFsDest = 6
        
      CASE (EL_Q2)

        ! Manually calculate the coordinates of the corners/midpoints on
        ! the reference element.
        Dxi(1,1)  =  -1.0_DP
        Dxi(1,2)  =  -1.0_DP

        Dxi(2,1)  =  1.0_DP
        Dxi(2,2)  =  -1.0_DP
        
        Dxi(3,1)  =  1.0_DP
        Dxi(3,2)  =  1.0_DP
        
        Dxi(4,1)  =  -1.0_DP
        Dxi(4,2)  =  1.0_DP

        Dxi(5,1)  =  0.0_DP
        Dxi(5,2)  =  -1.0_DP

        Dxi(6,1)  =  1.0_DP
        Dxi(6,2)  =  0.0_DP

        Dxi(7,1)  =  0.0_DP
        Dxi(7,2)  =  1.0_DP

        Dxi(8,1)  =  -1.0_DP
        Dxi(8,2)  =  0.0_DP
        
        Dxi(9,1)  =  0.0_DP
        Dxi(9,2)  =  0.0_DP
        
        nlocalDOFsDest = 9
        
      CASE DEFAULT
        CALL output_line ('Unsupported FE space in destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGrad2DInterpP1Q1cnf')
        CALL sys_halt()
      END SELECT

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      j = elem_igetCoordSystem(p_elementDistribution%itrialElement)
      
      ! Allocate memory and get local references to it.
      ! We abuse the system of cubature points here for the evaluation.
      CALL domint_initIntegration (rintSubset,nelementsPerBlock,nlocalDOFsDest,j,&
          p_rtriangulation%ndim,NVE)
      p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
      p_DcubPtsReal => rintSubset%p_DcubPtsReal
      p_Djac =>        rintSubset%p_Djac
      p_Ddetj =>       rintSubset%p_Ddetj
      p_Dcoords =>     rintSubset%p_DCoords

      ! Destination space is either P1, Q1, P2 or Q2. 
      ! Put the cubature point coordinates in the right format to the
      ! cubature-point array.
      ! Initialise all entries in p_DcubPtsRef with the same coordinates -
      ! as the cubature point coordinates are identical on all elements
      DO j=1,SIZE(p_DcubPtsRef,3)
        DO i=1,nlocalDOFsDest
          DO k=1,SIZE(p_DcubPtsRef,1)
            ! Could be solved using the TRANSPOSE operator - but often is's 
            ! faster this way...
            p_DcubPtsRef(k,i,j) = Dxi(i,k)
          END DO
        END DO
      END DO
      
      ! Allocate memory for the DOF's of all the elements.
      ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
      ALLOCATE(IdofsDest(indofDest,nelementsPerBlock))

      ! Allocate memory for the values of the derivatives in the corners
      ALLOCATE(Dderivatives(nlocalDOFsDest,nelementsPerBlock,2))

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
                     
      ! Loop over the elements - blockwise.
      DO IELset = 1, p_rtriangulation%NEL, PPGRD_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(p_rtriangulation%NEL,IELset-1+PPGRD_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        CALL dof_locGlobMapping_mult(p_rdiscrSource, p_IelementList(IELset:IELmax), &
                                     .TRUE.,IdofsTrial)

        ! Also calculate the global DOF's in our destination vector(s)
        CALL dof_locGlobMapping_mult(p_rdiscrDest, p_IelementList(IELset:IELmax), &
                                     .TRUE.,IdofsDest)

        ! We have the coordinates of the cubature points saved in the
        ! coordinate array from above. Unfortunately for nonparametric
        ! elements, we need the real coordinate.
        ! Furthermore, we anyway need the coordinates of the element
        ! corners and the Jacobian determinants corresponding to
        ! all the points.
        !
        ! At first, get the coordinates of the corners of all the
        ! elements in the current set. 
        
        CALL trafo_getCoords_sim (elem_igetTrafoType(&
            p_elementDistribution%itrialElement),&
            p_rtriangulation,p_IelementList(IELset:IELmax),p_Dcoords)
        
        ! Depending on the type of transformation, we must now choose
        ! the mapping between the reference and the real element.
        ! In case we use a nonparametric element, we need the 
        ! coordinates of the points on the real element, too.
        CALL trafo_calctrafo_sim (&
              p_rdiscrSource%RelementDistribution(icurrentElementDistr)%ctrafoType,&
              IELmax-IELset+1,nlocalDOFsDest,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1),&
              p_DcubPtsReal)
      
        ! Prepare the call to the evaluation routine of the analytic function.    
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
    
        ! At this point, we calculate the gradient information.
        ! As this routine handles the 2D case, we have an X- and Y-derivative.
        ! Calculate the X-derivative in the corners of the elements
        ! into Dderivatives(:,:,1) and the Y-derivative into
        ! Dderivatives(:,:,2).
        
        CALL fevl_evaluate_sim (rvector, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_elementDistribution%itrialElement, IdofsTrial, &
              nlocalDOFsDest, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV_X,&
              Dderivatives(:,1:IELmax-IELset+1_I32,1))        

        CALL fevl_evaluate_sim (rvector, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_elementDistribution%itrialElement, IdofsTrial, &
              nlocalDOFsDest, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV_Y,&
              Dderivatives(:,1:IELmax-IELset+1_I32,2))        
                    
         ! Sum up the derivative values in the destination vector.
         ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
         ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
         ! local DOF's in the destination space -- in that order!
                    
         DO i=1,IELmax-IELset+1
           DO j=1,nlocalDOFsDest
             p_DxDeriv(IdofsDest(j,i)) = p_DxDeriv(IdofsDest(j,i)) + Dderivatives(j,i,1)
             p_DyDeriv(IdofsDest(j,i)) = p_DyDeriv(IdofsDest(j,i)) + Dderivatives(j,i,2)
             
             ! Count how often a DOF was touched.
             p_IcontributionsAtDOF(IdofsDest(j,i)) = &
                p_IcontributionsAtDOF(IdofsDest(j,i))+1
           END DO
         END DO
                    
      END DO ! IELset
      
      ! Release memory
      CALL domint_doneIntegration(rintSubset)

      DEALLOCATE(Dderivatives)
      DEALLOCATE(IdofsDest)
      DEALLOCATE(IdofsTrial)

    END DO ! icurrentElementDistr

    ! We are nearly done. The final thing: divide the calculated derivatives by the
    ! number of elements adjacent to each vertex. That closes the calculation
    ! of the 'mean' of the derivatives.
    DO i=1,SIZE(p_DxDeriv)
      ! Div/0 should not occur, otherwise the triangulation is crap as there's a point
      ! not connected to any element!
      p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
      p_DyDeriv(i) = p_DyDeriv(i) / p_IcontributionsAtDOF(i)
    END DO
    
    ! Release temp data
    CALL storage_free (h_IcontributionsAtDOF)

  END SUBROUTINE

END MODULE
