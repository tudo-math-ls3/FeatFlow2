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
!#     -> Standard implementation for the calculation of the recovered 
!#        gradient of a scalar finite element function.
!#        A parameter admits to choose a method which is used with
!#        default parameters.
!#
!# Auxiliary routines, called internally.
!#
!# 1.) ppgrd_calcGrad2DInterpP12Q12cnf
!#     -> Calculate the reconstructed gradient as P1, Q1, P2 or Q2
!#        vector for an arbitrary conformal discretisation, 2D.
!#        Uses 1st order interpolation to reconstruct the gradient.
!#
!# 2.) ppgrd_calcGradSuperPatchRecov
!#     -> Calculate the reconstructed gradient as P1, Q1, P2 or Q2
!#        vector for an arbitrary conformal discretisation.
!#        Uses the superconvergent patch recovery technique suggested
!#        by Zienkiewicz and Zhu. 1D, 2D, 3D.
!#
!# 3.) ppgrd_calcGradLimAvgP1Q1cnf
!#     -> Calculate the reconstructed gradient as ~P1 or ~Q1 vector
!#        for an arbitrary conformal discretisation.
!#        Uses the limited gradient averaging technique by M. Möller
!#        and D. Kuzmin.
!#
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

!<constantblock description = "Identifiers for the type of patch used to recover the gradient vector.">

  ! Node-based patch: Use elements surrounding a particular node
  INTEGER, PARAMETER :: PPGRD_NODEPATCH = 0

  ! Element-based patch: Use elements surrounding a particular element
  INTEGER, PARAMETER :: PPGRD_ELEMPATCH = 1

  ! Face-based patch: Use subset of element-based patch which has common face
  INTEGER, PARAMETER :: PPGRD_FACEPATCH = 2

!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! Number of elements to handle simultaneously when building vectors
  INTEGER :: PPGRD_NELEMSIM   = 1000

  ! Number of patches to handle simultaneously when performing gradient recovery
  INTEGER :: PPGRD_NPATCHSIM  = 100
  
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
      ! Currently, the destination vector must be either pure Q1 or pure P1 or
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
        ! 1st order gradient
        CALL ppgrd_calcGrad2DInterpP12Q12cnf (rvectorScalar,rvectorGradient)
      CASE (PPGRD_ZZTECHNIQUE)
        ! 2nd order gradient with ZZ.
        ! Standard method is 'nodewise'.
        CALL ppgrd_calcGradSuperPatchRecov (rvectorScalar,rvectorGradient,PPGRD_NODEPATCH)
      END SELECT

    CASE DEFAULT
      ! Use ZZ as default.
      CALL ppgrd_calcGradSuperPatchRecov (rvectorScalar,rvectorGradient,PPGRD_NODEPATCH)
      !CALL output_line ('Unsupported dimension!',&
      !    OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      !CALL sys_halt()
    END SELECT

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppgrd_calcGrad2DInterpP12Q12cnf (rvector,rvectorGradient)

!<description>
  ! Calculates the recovered gradient of a scalar finite element function
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
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: IdofsTrial
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: IdofsDest
    
    ! Pointers to the X- and Y-derivative vector
    REAL(DP), DIMENSION(:), POINTER :: p_DxDeriv, p_DyDeriv
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

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
    
      ! If the element distribution is empty, skip it
      IF (p_elementDistribution%NEL .EQ. 0) CYCLE
    
      ! Get the number of local DOF's for trial functions
      ! in the source and destination vector.
      indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
      indofDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_elementDistribution%itrialElement)
      
      ! Initialise the cubature formula,
      ! That's a special trick here! The FE space of the destination vector
      ! is either P1 or Q1. We create the gradients in the corners of these points
      ! by taking the mean of the gradients of the source vector!
      !
      ! For this purpose, we initialise the 'trapezoidal rule' as cubature
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

      ! Get the number of elements in the element distribution.
      NEL = p_elementDistribution%NEL
      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPGRD_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPGRD_NELEMSIM)
      
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

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppgrd_calcGradSuperPatchRecov (rvectorScalar,rvectorGradient,cpatchType)

!<description>
    ! Calculates the recovered gradient of a scalar finite element function
    ! by means of the superconvergent patch recovery technique suggested
    ! by Zienkiewicz and Zhu. Supports conformal discretisations in arbitrary
    ! spatial dimensions with $P_1$, $Q_1$, $P_2$ and $Q_2$ finite elements
    ! mixed in the source and destination vectors.
!</description>

!<input>
    ! The FE solution vector. Represents a scalar FE function.
    TYPE(t_vectorScalar), INTENT(IN)         :: rvectorScalar
    
    ! The type of patch used to recover the gradient values
    INTEGER, INTENT(IN)                      :: cpatchType
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
    INTEGER :: IVE,NVE,NVEMax
    INTEGER :: icurrentElementDistr,ilastElementDistr,ilocalElementDistr
    INTEGER :: i,j,k,ipoint,idx,icoordSystem
    INTEGER :: IPATCH,NPATCH,PATCHset,PATCHmax
    LOGICAL :: bnonparTrial
    INTEGER(PREC_ELEMENTIDX)    :: IEL,JEL,KEL
    INTEGER(PREC_VERTEXIDX)     :: IVT
    REAL(DP), DIMENSION(NDIM3D) :: Dval

    ! Array to tell the element which derivatives to calculate
    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder,BderDest

    ! Cubature point coordinates on the reference element
    REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega

    ! Flag if discretisation is uniform
    LOGICAL :: bisuniform
    
    ! Number of cubature points on the reference element
    INTEGER :: ncubp

    ! Maximum number of cubature points on the reference element in the block
    INTEGER :: ncubpMax
    
    ! Number of local degees of freedom for test functions
    INTEGER :: indofTrial

    ! Maximum numer of local degrees of freedom for test functiions in the block
    INTEGER :: indofTrialMax

    ! Number of local degrees of freedom for destination space
    INTEGER :: indofDest

    ! Maximum number of local degrees of freedom for destination space
    INTEGER :: indofDestMax

    ! Number of 'cubature points'. As we acutally don't do cubature
    ! here, this coincides with the number of DOF's on each element
    ! in the destination space.
    INTEGER :: nlocalDOFsDest

    ! Maximum number of 'cubature points'. As we acutally don't do 
    ! cubature here, this coincides with the number of DOF's on each 
    ! element in the destination space.
    INTEGER :: nlocalDOFsDestMax

    ! Total number of sampling/interpolation points in set of patches
    INTEGER :: nspoints,nipoints

    ! The triangulation structure - to shorten some things...
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! The spatial discretisation structure - to shorten some things...
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrSource
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrDest

    ! Current element distribution in use
    TYPE(t_elementDistribution), POINTER :: p_elementDistribution
    TYPE(t_elementDistribution), POINTER :: p_elementDistrDest

    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    TYPE(t_domainIntSubset) :: rintSubset
    TYPE(t_domainIntSubset) :: rintSubsetDest

    ! An allocatable array accepting the starting positions of each patch
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), ALLOCATABLE :: IelementsInPatchIdx

    ! An allocatable array accepting the element numbers of each patch
    ! Note that the first member for each patch defines the patch details, i.e.,
    ! the cubature formular that should be used for this patch, etc.
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), ALLOCATABLE :: IelementsInPatch

    ! An allocatable array accepting the number of corners per element
    INTEGER, DIMENSION(:), ALLOCATABLE :: IelementNVEInPatch

    ! An allocatable array accepting the number of cubature points per element
    INTEGER, DIMENSION(:), ALLOCATABLE :: IelementNcubpInPatch

    ! An allocatable array accepting the number of sampling points per patch
    INTEGER, DIMENSION(:), ALLOCATABLE :: Inpoints
    
    ! An allocatable array accepting the DOF's of a set of elements.
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: IdofsTrial
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: IdofsDest

    ! An allocatable array accepting the coordinates of the patch bounding group
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: DpatchBound

    ! An allocatable array accepting the polynomials of the set of elements.
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Dpolynomials
    REAL(DP), DIMENSION(:),       ALLOCATABLE :: DpolynomialsMixed

    ! An allocatable array accepting the values of the FE functions
    ! that are computed by the callback routine.
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: DcoefficientsMixed

    ! An allocatable array accepting the averaged gradient values
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dderivatives

    ! Pointer to element-at-vertex index array of the triangulation
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx

    ! Pointer to element-at-vertex list of the triangulation
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex

    ! Pointer to elements-at-element list of the triangulation
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement

     ! Pointer to vertices-at-element list of the triangulation
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement

    ! Pointer to an element-number list
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList

    ! Pointer to vertex coordinates of the triangulation
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords

    ! Pointer to element distribution identifier list.
    INTEGER, DIMENSION(:), POINTER :: p_IelementDistr

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

    ! Number of patches in a block
    INTEGER :: npatchesPerBlock

    ! Number of patches currently blocked
    INTEGER :: npatchesInCurrentBlock

    ! Number of elements in a block
    INTEGER :: nelementsPerBlock

    ! Pointer to an array that counts the number of elements adjacent to a vertex.
    ! Ok, there's the same information in the triangulation, but that's not
    ! based on DOF's! Actually, we'll calculate how often we touch each DOF 
    ! in the destination space.
    INTEGER :: h_IcontributionsAtDOF
    INTEGER(I32), DIMENSION(:), POINTER :: p_IcontributionsAtDOF
    
    ! Pointers to the X-, Y- and Z-derivative vector
    REAL(DP), DIMENSION(:), POINTER :: p_DxDeriv, p_DyDeriv, p_DzDeriv

    ! Auxiliary integers
    INTEGER :: icubp,idim
    
    !-----------------------------------------------------------------------
    ! (0)  Initialisation
    !-----------------------------------------------------------------------
       
    ! Get the discretisation structures of the source- and destination space.
    ! Note that we assume here that all derivatives are discretised the same way!
    p_rdiscrSource => rvectorScalar%p_rspatialDiscretisation
    p_rdiscrDest   => rvectorGradient%p_rblockDiscretisation%RspatialDiscretisation(1)
    
    ! Check if we have a non-uniform discretisation structure and if nodal 
    ! patches should be used. This does not make too much sense since it is
    ! not clear which type of element should be adopted for the patch elements.
    !
    ! Theoretically, one could allow for an optional parameter which defines
    ! the element type of the "patch" elements but this requires more complicated
    ! code. There are other possibilities for the patches that should be used !!!
    IF ((p_rdiscrSource%ccomplexity .NE. SPDISC_UNIFORM .OR. &
         p_rdiscrDest%ccomplexity   .NE. SPDISC_UNIFORM) .AND. &
         cpatchType .EQ. PPGRD_NODEPATCH) THEN 
      CALL output_line('Nodal patches are not available for non-uniform discretisations!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      CALL sys_halt()
    END IF

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscrSource%p_rtriangulation
    
    ! Get a pointer to the vertices-at-element and vertex coordinates array
    CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    CALL storage_getbase_double2D (p_rtriangulation%h_DvertexCoords,  p_DvertexCoords)

    ! Get pointers to the derivative destination vector.
    SELECT CASE(p_rtriangulation%ndim)
    CASE (NDIM1D)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      CALL lalg_clearVectorDble (p_DxDeriv)

    CASE (NDIM2D)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      CALL lalg_clearVectorDble (p_DxDeriv)
      CALL lalg_clearVectorDble (p_DyDeriv)

    CASE (NDIM3D)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(3),p_DzDeriv)
      CALL lalg_clearVectorDble (p_DxDeriv)
      CALL lalg_clearVectorDble (p_DyDeriv)
      CALL lalg_clearVectorDble (p_DzDeriv)

    CASE DEFAULT
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      CALL sys_halt()
    END SELECT
      
    ! Array that allows the calculation about the number of elements
    ! meeting in a vertex, based onm DOF's.
    CALL storage_new ('ppgrd_calcGradSuperPatchRecov','DOFContrAtVertex',&
        dof_igetNDofGlob(p_rdiscrDest), ST_INT, h_IcontributionsAtDOF, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(h_IcontributionsAtDOF, p_IcontributionsAtDOF)
        
    ! We only need the derivatives of trial functions
    Bder = .FALSE.
    SELECT CASE (p_rtriangulation%ndim)
    CASE (NDIM1D)
      Bder(DER_DERIV1D_X) = .TRUE.

    CASE (NDIM2D)
      Bder(DER_DERIV2D_X) = .TRUE.
      Bder(DER_DERIV2D_Y) = .TRUE.

    CASE (NDIM3D)
      Bder(DER_DERIV3D_X) = .TRUE.
      Bder(DER_DERIV3D_Y) = .TRUE.
      Bder(DER_DERIV3D_Z) = .TRUE.

    CASE DEFAULT
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      CALL sys_halt()
    END SELECT
    
    ! For the recovery, we only need the function values of trial functions
    BderDest = .FALSE.
    BderDest(DER_FUNC) = .TRUE.
    

    ! Do we have a uniform triangulation? Would simplify a lot...
    bisUniform = (p_rdiscrSource%ccomplexity .EQ. SPDISC_UNIFORM) .AND. &
                 (p_rdiscrDest%ccomplexity   .EQ. SPDISC_UNIFORM)

    IF (.NOT. bisUniform) THEN
      ! Things are more complicated if one of the discretisations is not uniform.
      ! In this case, we always have to consider the maximum number of quadrature
      ! points, the largest number of local DOF's, etc. 

      indofTrialMax     = 0
      indofDestMax      = 0
      NVEmax            = 0
      ncubpMax          = 0
      nlocalDOFsDestMax = 0
      
      ! Set pointer to list of element distributions which exists in this case
      CALL storage_getbase_int(p_rdiscrSource%h_IelementDistr, p_IelementDistr)

      ! Loop over all element distributions of the source discretisation
      DO icurrentElementDistr = 1, p_rdiscrSource%inumFESpaces
        
        ! Activate the current element distribution
        p_elementDistribution => p_rdiscrSource%RelementDistribution(icurrentElementDistr)
        
        ! Cancel if this element distribution is empty.
        IF (p_elementDistribution%NEL .EQ. 0) CYCLE
        
        ! Get the number of local DOF's for trial functions
        indofTrial    = elem_igetNDofLoc(p_elementDistribution%itrialElement)
        indofTrialMax = MAX(indofTrialMax,indofTrial)
        
        ! Get the number of corner vertices of the element
        NVE    = elem_igetNVE(p_elementDistribution%itrialElement)
        NVEmax = MAX (NVEmax,NVE)
        
        ! Get cubature weights and point coordinates on the reference element
        CALL cub_getCubPoints(p_elementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
        ncubpMax = MAX(ncubpMax,ncubp)
      END DO

      ! Loop over all element distributions of the destination discretisation
      DO icurrentElementDistr = 1, p_rdiscrDest%inumFESpaces
        
        ! Activate the current element distribution
        p_elementDistribution => p_rdiscrDest%RelementDistribution(icurrentElementDistr)

        ! Cancel if this element distribution is empty.
        IF (p_elementDistribution%NEL .EQ. 0) CYCLE

        ! Get the number of local DOF's for trial functions
        indofDest    = elem_igetNDofLoc(p_elementDistribution%itrialElement)
        indofDestMax = MAX(indofDestMax,indofDest)

        ! Get cubature weights and point coordinates on the reference element
        CALL cub_getCubPoints(p_elementDistribution%ccubTypeEval, nlocalDOFsDest, Dxi, Domega)
        nlocalDOFsDestMax = MAX(nlocalDOFsDestMax,nlocalDOFsDest)
      END DO
    END IF
        

    ! Set pointers which are used to assemble patches
    SELECT CASE(cpatchType)
    CASE (PPGRD_NODEPATCH)
      ! Get the elements-at-vertex index array.
      CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
      ! Get the elements-at-vertex array.
      CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
            
    CASE (PPGRD_ELEMPATCH)
      ! Get the elements-at-vertex index array.
      CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
      ! Get the elements-at-vertex array.
      CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
      ! Get the elements-at-element array.
      CALL storage_getbase_int2D (p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
      ! Get the vertices-at-element array.
      CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)     

    CASE (PPGRD_FACEPATCH)
      ! Get the elements-at-element array.
      CALL storage_getbase_int2D (p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)

    CASE DEFAULT
      CALL output_line('Invalid patch type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGrad2DSuperPatchRecov')
      CALL sys_halt()
    END SELECT


    !---------------------------------------------------------------------------
    ! (1)  Superconvergent patch recovery - main loop
    !---------------------------------------------------------------------------
    
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation. Note that
    ! nodal patches are only allowed for uniform discretisations, that is,
    ! we can be sure, that there is exactly one element distribution if
    ! nodal patches should be considered.

    DO icurrentElementDistr = 1, p_rdiscrSource%inumFESpaces

      ! Activate the current element distribution
      p_elementDistribution => p_rdiscrSource%RelementDistribution(icurrentElementDistr)
      
      ! If the element distribution is empty, skip it
      IF (p_elementDistribution%NEL .EQ. 0) CYCLE
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_elementDistribution%itrialElement)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_elementDistribution%h_IelementList, p_IelementList)
      
      ! Get number of patches in current element distribution
      SELECT CASE(cpatchType)
      CASE (PPGRD_NODEPATCH)
        ! Recall that for nodal-baes patches there MUST be exactly one element 
        ! distribution so that the number of patches equals the total number
        ! of vertices in the whole triangulation !!!
        NPATCH = p_rtriangulation%NVT
        
      CASE (PPGRD_ELEMPATCH,PPGRD_FACEPATCH)
        ! The number of patches equals the number of elements 
        ! in the current element distribution.
        NPATCH = p_elementDistribution%NEL
        
      CASE DEFAULT
        CALL output_line('Invalid patch type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGrad2DSuperPatchRecov')
        CALL sys_halt()
      END SELECT

      ! Determine actual number of patches in this block
      npatchesPerBlock = MIN(PPGRD_NPATCHSIM,NPATCH)

      ! Allocate memory for element numbers in patch index array
      ALLOCATE(IelementsInPatchIdx(npatchesPerBlock+1))
      
      ! Allocate memory for coordinates of patch bounding groups
      ALLOCATE(DpatchBound(p_rtriangulation%ndim,2,npatchesPerBlock))

      ! Allocate memory for number of sampling points
      IF (.NOT. bisUniform) ALLOCATE(Inpoints(npatchesPerBlock))
      
      ! Loop over the patches - blockwise.
      DO PATCHset = 1, NPATCH , PPGRD_NPATCHSIM
        
        !-------------------------------------------------------------------------
        ! Phase 1: Determine all elements in the set of patches
        !-------------------------------------------------------------------------
        
        ! We always handle PPGRD_NPATCHSIM patches simultaneously.
        ! How many patches have we actually here?
        ! Get the maximum patch number, such that we handle at most 
        ! PPGRD_NPATCHSIM patches simultaneously.
        PATCHmax = MIN(NPATCH,PATCHset-1+PPGRD_NPATCHSIM)
        
        ! Calculate the number of patches currently blocked
        npatchesInCurrentblock = PATCHmax-PATCHset+1
        
        ! Depending on the patch type, the element numbers that are contained in
        ! each patch are different and must be determined from different structures
        SELECT CASE(cpatchType)
        CASE (PPGRD_NODEPATCH)
          
          ! We actually know how many elements are adjacent to each node, so we can
          ! allocate memory for element numbers in patch without 'dummy' loop
          nelementsPerBlock = p_IelementsAtVertexIdx(PATCHmax+1)-p_IelementsAtVertexIdx(PATCHset)
          ALLOCATE(IelementsInPatch(nelementsPerBlock))
          
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Loop over the patches in the set
          DO IPATCH = 1, npatchesInCurrentBlock
            
            ! Store index of first element in this patch
            IelementsInPatchIdx(IPATCH) = nelementsPerBlock + 1
            
            ! Get global vertex number
            IVT = PATCHset-1+IPATCH
            
            ! Loop over elements adjacent to vertex
            DO idx = p_IelementsAtVertexIdx(IVT),p_IelementsAtVertexIdx(IVT+1)-1
              
              ! Get element number
              IEL = p_IelementsAtVertex(idx)
              
              ! Increase element counter
              nelementsPerBlock = nelementsPerBlock + 1
              IelementsInPatch(nelementsPerBlock) = IEL
            END DO
          END DO
          
          ! Store index of last element in last patch increased by one
          IelementsInPatchIdx(npatchesInCurrentBlock+1) = nelementsPerBlock + 1
         
          
        CASE (PPGRD_ELEMPATCH)
          
          ! Unfortunately, we do not directly know how many elements are in the neighbourhood
          ! of each element. But we know, how many elements are adjacent to each corner of
          ! a particular element. If we sum up these numbers we get an upper bound for the
          ! number of elements present in each patch. Obviously, the center element is multiply
          ! counted. Moreover, all edge/face neighbours are also counted twice.
          
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Ok, let's do a dummy loop to determine the number of elements in the block
          DO IPATCH = 1, npatchesInCurrentBlock
            
            ! Get the global element number from the list of elements in distribution
            IEL = p_IelementList(PATCHset+IPATCH-1)
            
            ! Loop over corner nodes
            DO ive =1, NVE
              
              ! Get global vertex number
              IVT = p_IverticesAtElement(ive,IEL)
              
              ! Count number of elements surrounding corner node
              nelementsPerBlock = nelementsPerBlock + &
                  (p_IelementsAtVertexIdx(IVT+1)-p_IelementsAtVertexIdx(IVT))
            END DO
            
            ! Ok, we counted element IEL NVE-times but it is only required once
            nelementsPerBlock = nelementsPerBlock - (NVE-1)
            
            ! Moreover, each adjacent element is counted twice if it is not the boundary
            nelementsPerBlock = nelementsPerBlock - COUNT(p_IneighboursAtElement(1:NVE,IEL) > 0)
          END DO
          
          ! That's it, we can allocate memory for elements numbers
          ALLOCATE(IelementsInPatch(nelementsPerBlock))
          
          ! Now, we have to fill it with the element numbers
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Loop over the patches in the set
          DO IPATCH = 1, npatchesInCurrentBlock
            
            ! Store index of first element in this patch
            IelementsInPatchIdx(IPATCH) = nelementsPerBlock + 1
            
            ! Get the global element number from the list of elements in distribution
            IEL = p_IelementList(PATCHset+IPATCH-1)
            
            ! Do not forget to store the element itself
            nelementsPerBlock = nelementsPerBlock + 1
            IelementsInPatch(nelementsPerBlock) = IEL
            
            ! Loop over corner nodes
            DO ive = 1, NVE
              
              ! Get the globale vertex number
              IVT = p_IverticesAtElement(ive,IEL)
              
              ! Get element number adjacent to element IEL
              JEL = p_IneighboursAtElement(ive,IEL)
              
              ! Loop over elements adjacent to vertex
              DO idx = p_IelementsAtVertexIdx(IVT),p_IelementsAtVertexIdx(IVT+1)-1
                
                ! Get element number
                KEL = p_IelementsAtVertex(idx)
                
                ! Do not consider KEL = IEL and KEL = JEL to prevent the same
                ! element number to be present in a patch multiple times
                IF (KEL .EQ. IEL .OR. KEL .EQ. JEL) CYCLE
                
                ! Increase element counter
                nelementsPerBlock = nelementsPerBlock + 1
                IelementsInPatch(nelementsPerBlock) = KEL
              END DO
            END DO
          END DO
          
          ! Store index of last element in last patch increased by one
          IelementsInPatchIdx(npatchesInCurrentBlock+1) = nelementsPerBlock + 1
                    
        CASE (PPGRD_FACEPATCH)
          
          ! We actually know how many elements are adjacent to each element, so we
          ! can allocate memory for element numbers in patch without 'dummy' loop
          nelementsPerBlock = (NVE+1)*npatchesPerBlock
          ALLOCATE(IelementsInPatch(nelementsPerBlock))
          
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Loop over the patches in the set
          DO IPATCH = 1, npatchesInCurrentBlock
            
            ! Store index of first element in this patch
            IelementsInPatchIdx(IPATCH) = nelementsPerBlock + 1
            
            ! Get the global element number from the list of elements in distribution
            IEL = p_IelementList(PATCHset+IPATCH-1)
            
            ! Do not forget to store the element itself
            nelementsPerBlock = nelementsPerBlock + 1
            IelementsInPatch(nelementsPerBlock) = IEL
            
            ! Loop over adjacent elements
            DO ive = 1, NVE
              
              ! Get element number adjacent to element IEL
              JEL = p_IneighboursAtElement(ive,IEL)
              
              ! Check if element neighbour is the boundary, then skip it
              IF (JEL .EQ. 0) CYCLE
              
              ! Increase element counter
              nelementsPerBlock = nelementsPerBlock + 1
              IelementsInPatch(nelementsPerBlock) = JEL
            END DO
          END DO
          
          ! Store index of last element in last patch increased by one
          IelementsInPatchIdx(npatchesInCurrentBlock+1) = nelementsPerBlock + 1
          
        CASE DEFAULT
          CALL output_line('Invalid patch type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGrad2DSuperPatchRecov')
          CALL sys_halt()
        END SELECT
        
        
        !-----------------------------------------------------------------------
        ! Phase 2: Perform Least-squares fitting on the set of patches
        !-----------------------------------------------------------------------
        
        ! Do we have a uniform discretisation? Would simplify a lot...
        IF (bisUniform) THEN
          
          ! Yes, the discretisation is uniform. In this case, we have to set the 
          ! pointers, dimensions, etc. just once and can work blockwise.
          
          ! Active element distribution
          p_elementDistribution => p_rdiscrSource%RelementDistribution(icurrentElementDistr)
          p_elementDistrDest    => p_rdiscrDest%RelementDistribution(icurrentElementDistr)
          
          ! Get the number of local DOF's for trial functions
          indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
          indofDest  = elem_igetNDofLoc(p_elementDistrDest%itrialElement)
          
          ! Get the number of corner vertices of the element
          NVE = elem_igetNVE(p_elementDistribution%itrialElement)
          

          !---------------------------------------------------------------------
          ! Step 1:  Prepare the source FE space
          !---------------------------------------------------------------------
          
          ! Allocate memory for the DOF's of all the elements
          ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
          
          ! Calculate the global DOF's into IdofsTrial.
          CALL dof_locGlobMapping_mult(p_rdiscrSource, &
              IelementsInPatch(1:nelementsPerBlock), .FALSE., IdofsTrial)
          
          ! Initialise the cubature formula for the source element distribution
          ! Get cubature weights and point coordinates on the reference element
          CALL cub_getCubPoints(p_elementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
          
          ! Get from the trial element space the type of coordinate system
          ! that is used there:
          icoordSystem = elem_igetCoordSystem(p_elementDistribution%itrialElement)
          
          ! Allocate memory and get local references to it. This domain integration 
          ! structure stores all information of the source FE space. 
          CALL domint_initIntegration (rintSubset, nelementsPerBlock, &
              ncubp, icoordSystem, p_rtriangulation%ndim, NVE)
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
          
          ! Check if one of the trial/test elements is nonparametric
          bnonparTrial = elem_isNonparametric(p_elementDistribution%itestElement)
          
          ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
          ! p_DcubPtsRef - depending on whether the space is parametric or not.
          IF (bnonparTrial) THEN
            p_DcubPtsTrial => p_DcubPtsReal
          ELSE
            p_DcubPtsTrial => p_DcubPtsRef
          END IF
          
          ! We have the coordinates of the cubature points saved in the
          ! coordinate array from above. Unfortunately for nonparametric
          ! elements, we need the real coordinate.
          ! Furthermore, we anyway need the coordinates of the element
          ! corners and the Jacobian determinants corresponding to
          ! all the points.
          
          ! At first, get the coordinates of the corners of all the
          ! elements in the current set of elements.
          CALL trafo_getCoords_sim (&
              elem_igetTrafoType(p_elementDistribution%itrialElement), &
              p_rtriangulation, IelementsInPatch(1:nelementsPerBlock), p_Dcoords)
          
          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType, &
              nelementsPerBlock, ncubp, p_Dcoords, p_DcubPtsRef, p_Djac, &
              p_Ddetj, p_DcubPtsReal)
          
          !---------------------------------------------------------------------
          ! Step 2:  Perform sampling of consistent gradient values
          !---------------------------------------------------------------------
          
          ! Allocate memory for the values of the derivaties in the corners
          ALLOCATE(Dcoefficients(ncubp,nelementsPerBlock,p_rtriangulation%ndim))
          
          ! Calculate the derivative of the FE function in the cubature
          ! points: u_h(x,y,z) and save the result to Dcoefficients(:,:,1..3)
          SELECT CASE(p_rtriangulation%ndim)
          CASE (NDIM1D)
            CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_elementDistribution%itrialElement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV1D_X, Dcoefficients(:,:,1))
            
          CASE (NDIM2D)
            CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_elementDistribution%itrialElement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV2D_X, Dcoefficients(:,:,1))
            
            CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_elementDistribution%itrialElement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV2D_Y, Dcoefficients(:,:,2))
            
          CASE (NDIM3D)
            CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_elementDistribution%itrialElement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV3D_X, Dcoefficients(:,:,1))
            
            CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_elementDistribution%itrialElement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV3D_Y, Dcoefficients(:,:,2))
            
            CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_elementDistribution%itrialElement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV3D_Z, Dcoefficients(:,:,3))

          CASE DEFAULT
            CALL output_line('Invalid spatial dimension!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
            CALL sys_halt()
          END SELECT
          

          !---------------------------------------------------------------------
          ! Step 3: Prepare least-squares fitting
          !---------------------------------------------------------------------
          
          ! First, get the coordinates of the bounding group for the set of
          ! elements present in each patch. In addition, calculate the coordinates
          ! of the corners of the constant Jacobian patch "elements".
          CALL calc_patchBoundingGroup_sim(IelementsInPatchIdx, npatchesInCurrentBlock, &
              NVE, p_Dcoords, DpatchBound(:,:,1:npatchesInCurrentBlock))
          
          ! Next, we need to convert the physical coordinates of the curvature points
          ! to the local coordinates of the constant Jacobian "patch" elements
          CALL calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, &
              npatchesInCurrentBlock, NVE, p_DcubPtsReal, DpatchBound, p_DcubPtsRef)
          
          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType, &
              nelementsPerBlock, ncubp, p_Dcoords, p_DcubPtsRef, p_Djac, p_Ddetj)
          
          ! Allocate memory for the patch interpolants matrices
          ! Note that the second dimension is DER_FUNC=1 and could be omitted. However,
          ! all routines for simultaneous element evaluation require a 4d-array so that
          ! the array Dpolynomials is artificially created as 4d-array.
          ALLOCATE(Dpolynomials(indofTrial,DER_FUNC,ncubp,nelementsPerBlock))
          
          ! Evaluate the trial functions of the constant Jacobian patch "element" for all
          ! cubature points of the elements present in the patch and store each polynomial 
          ! interpolation in the rectangular patch matrix used for least-squares fitting.
          CALL elem_generic_sim(p_elementDistribution%itrialElement,&
              p_Dcoords, p_Djac, p_Ddetj, BderDest, Dpolynomials, ncubp,&
              nelementsPerBlock, p_DcubPtsTrial)
          
          
          !-----------------------------------------------------------------------
          ! Step 4: Perform least-squares fitting
          !-----------------------------------------------------------------------
          
          ! Allocate memory for the derivative values
          ALLOCATE(Dderivatives(indofTrial,npatchesInCurrentBlock,p_rtriangulation%ndim))
          
          ! Compute the patch averages by solving $(P^T * P) * x = (P^T) * b$ for x
          CALL calc_patchAverages_sim(IelementsInPatchIdx, npatchesInCurrentBlock, ncubp, &
              indofTrial,  Dcoefficients, Dpolynomials,Dderivatives)
          

          !-----------------------------------------------------------------------
          ! Step 5: Prepare the destination FE space
          !-----------------------------------------------------------------------
          
          ! Allocate memory for the DOF's of all the elements
          ALLOCATE(IdofsDest(indofDest,nelementsPerBlock))
          
          ! Also calculate the global DOF's in our destination vector(s)
          CALL dof_locGlobMapping_mult(p_rdiscrDest, &
              IelementsInPatch(1:nelementsPerBlock), .FALSE., IdofsDest)
          
          ! Initialise the cubature formula. That's a special trick here!
          ! In particular, we only need the physical coordinates of the
          ! nodal evaluation points in the source vectors.
          ! For this purpose, we initialise the 'trapezoidal rule' as cubature
          ! formula. Ok, we are not using cubature at all here (thus ignoring
          ! the cubature weights completely), but that way we get the
          ! coordinates of the corners on the reference element automatically --
          ! which coincide with the points where we want to create the gradients!
          !
          ! Note: The returned nlocalDOFsDest will coincide with the number of local DOF's
          ! on each element indofDest!
          CALL calc_cubatureDest(&
              elem_getPrimaryElement(p_elementDistrDest%itrialElement), &
              nlocalDOFsDest, Dxi, Domega)
          
          ! Get from the trial element space the type of coordinate system
          ! that is used there:
          icoordSystem = elem_igetCoordSystem(p_elementDistrDest%itrialElement)
          
          ! Allocate memory and get local references to it. This domain integration 
          ! structure stores all information of the destination FE space.
          CALL domint_initIntegration (rintSubsetDest, nelementsPerBlock, &
              nlocalDOFsDest, icoordSystem, p_rtriangulation%ndim, NVE)
          p_DcubPtsRef =>  rintSubsetDest%p_DcubPtsRef
          p_DcubPtsReal => rintSubsetDest%p_DcubPtsReal
          p_Djac =>        rintSubsetDest%p_Djac
          p_Ddetj =>       rintSubsetDest%p_Ddetj
          p_Dcoords =>     rintSubsetDest%p_DCoords
          
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
          
          ! Check if one of the trial/test elements is nonparametric
          bnonparTrial = elem_isNonparametric(p_elementDistrDest%itestElement)
          
          ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
          ! p_DcubPtsRef - depending on whether the space is parametric or not.
          IF (bnonparTrial) THEN
            p_DcubPtsTrial => p_DcubPtsReal
          ELSE
            p_DcubPtsTrial => p_DcubPtsRef
          END IF
          
          ! We have the coordinates of the cubature points saved in the
          ! coordinate array from above. Unfortunately for nonparametric
          ! elements, we need the real coordinate.
          ! Furthermore, we anyway need the coordinates of the element
          ! corners and the Jacobian determinants corresponding to
          ! all the points.
          
          ! At first, get the coordinates of the corners of all the
          ! elements in the current set of elements.
          CALL trafo_getCoords_sim (elem_igetTrafoType(p_elementDistrDest%itrialElement), &
              p_rtriangulation, IelementsInPatch(1:nelementsPerBlock), p_Dcoords)
          
          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          CALL trafo_calctrafo_sim (p_elementDistrDest%ctrafoType, &
              nelementsPerBlock, nlocalDOFsDest, p_Dcoords, p_DcubPtsRef, p_Djac, &
              p_Ddetj, p_DcubPtsReal)
          
          ! Next, we need to convert the physical coordinates of the curvature points
          ! to the local coordinates of the constant Jacobian "patch" elements
          CALL calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, &
              npatchesInCurrentBlock, NVE,  p_DcubPtsReal, DpatchBound, p_DcubPtsRef)
          
          ! We do not need the corner coordinates of the elements in the destination
          ! FE space but that of the "patch" elements in the source FE space
          p_Dcoords => rintSubset%p_Dcoords

          ! Calculate the transformation from the reference elements to the real ones
          CALL trafo_calctrafo_sim (p_elementDistrDest%ctrafoType, &
              nelementsPerBlock, nlocalDOFsDest, p_Dcoords, p_DcubPtsRef, p_Djac, &
              p_Ddetj)
          
          ! Reallocate memory for the patch interpolants
          IF (ncubp .LT. nlocalDOFsDest) THEN
            DEALLOCATE(Dpolynomials)
            ALLOCATE(Dpolynomials(indofTrial,DER_FUNC,nlocalDOFsDest,nelementsPerBlock))
          END IF
          
          ! Evaluate the basis functions for the cubature points of the destination FE space
          CALL elem_generic_sim(p_elementDistribution%itrialElement,&
              p_Dcoords, p_Djac, p_Ddetj, BderDest, Dpolynomials, &
              nlocalDOFsDest, nelementsPerBlock, p_DcubPtsTrial)
          
          !---------------------------------------------------------------------
          ! Step 6: Evaluate the averaged derivative values at the cubature
          !         points of the destination FE space and scatter them to
          !         the global degrees of freedom
          !---------------------------------------------------------------------
          
          SELECT CASE (p_rtriangulation%ndim)
          CASE(NDIM1D) 
            ! Loop over the patches in the set
            DO ipatch = 1, npatchesInCurrentBlock

              ! Loop over elements in patch
              DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
                
                ! Loop over local degrees of freedom
                DO ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  DO j = 1, indofTrial
                    Dval(1) = Dval(1) + Dderivatives(j,ipatch,1) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  END DO
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                END DO
              END DO
            END DO

          CASE (NDIM2D)
            ! Loop over the patches in the set
            DO ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
                
                ! Loop over local degrees of freedom
                DO ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  DO j = 1, indofTrial
                    Dval(1:2) = Dval(1:2) + Dderivatives(j,ipatch,1:2) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  END DO
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                END DO
              END DO
            END DO

          CASE (NDIM3D)
            ! Loop over the patches in the set
            DO ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
                
                ! Loop over local degrees of freedom
                DO ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  DO j = 1, indofTrial
                    Dval(1:3) = Dval(1:3) + Dderivatives(j,ipatch,1:3) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  END DO
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  p_DzDeriv(IdofsDest(ipoint,idx)) = p_DzDeriv(IdofsDest(ipoint,idx)) + Dval(3)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                END DO
              END DO
            END DO
            
          CASE DEFAULT
            CALL output_line('Invalid spatial dimension!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
            CALL sys_halt()
          END SELECT
          
          ! Release memory
          CALL domint_doneIntegration(rintSubset)
          CALL domint_doneIntegration(rintSubsetDest)
          
          ! Deallocate temporary memory
          DEALLOCATE(IdofsDest)
          DEALLOCATE(Dderivatives)
          DEALLOCATE(Dpolynomials)
          DEALLOCATE(Dcoefficients)
          DEALLOCATE(IdofsTrial)
          
        ELSE
          
          ! No, the discretisation is not uniform. In this case, we have to check the type
          ! of element for each individual element and adopt the FE spaces accordingly.
          ilastElementDistr = 0

          ! Reset number of sampling points
          Inpoints = 0

          !-----------------------------------------------------------------------
          ! Step 0:  Allocate temporal memory; note that the dimensions of all
          !          arrays are set to the maximum values required in the process
          !-----------------------------------------------------------------------

          ! Allocate memory for number of corners per element
          ALLOCATE(IelementNVEInPatch(nelementsPerBlock))

          ! Allocate memory for number of cubature points per element
          ALLOCATE(IelementNcubpInPatch(nelementsPerBlock))
          
          ! Allocate memory for the DOF's of all the elements
          ALLOCATE(IdofsTrial(indofTrialMax,nelementsPerBlock))

          ! Allocate memory for the DOF's of all the elements
          ALLOCATE(IdofsDest(indofDestMax,nelementsPerBlock))

          ! Allocate memory for coefficient values; for mixed triangulations
          ! $ncubpMax$ is an upper bound for the number of cubature points.
          ALLOCATE(Dcoefficients(ncubpMax,nelementsPerBlock,p_rtriangulation%ndim))
!!$          ALLOCATE(DcoefficientsMixed(ncubpMax*nelementsPerBlock,p_rtriangulation%ndim))

          ! Calculate the global DOF's into IdofsTrial.
          CALL dof_locGlobMapping_mult(p_rdiscrSource, &
              IelementsInPatch(1:nelementsPerBlock), .FALSE., IdofsTrial)
          
          ! Also calculate the global DOF's in our destination vector(s)
          CALL dof_locGlobMapping_mult(p_rdiscrDest, &
              IelementsInPatch(1:nelementsPerBlock), .FALSE., IdofsDest)

          ! Allocate memory and get local references to it. This domain integration 
          ! structure stores all information of the source FE space. 
          CALL domint_initIntegration (rintSubset, nelementsPerBlock, &
              ncubpMax, TRAFO_CS_BARY2DTRI, p_rtriangulation%ndim, NVEMax)
          p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
          p_DcubPtsReal => rintSubset%p_DcubPtsReal
          p_Djac =>        rintSubset%p_Djac
          p_Ddetj =>       rintSubset%p_Ddetj
          p_Dcoords =>     rintSubset%p_DCoords

          ! Initialise arrays with zeros
          p_DcubPtsRef  = 0
          p_DcubPtsReal = 0
          p_Djac        = 0
          p_Ddetj       = 0
          p_Dcoords     = 0

          ! Allocate memory. This domain integration structure stores 
          ! all information of the destination FE space. 
          CALL domint_initIntegration (rintSubsetDest, nelementsPerBlock, &
              nlocalDOFsDestMax, TRAFO_CS_BARY2DTRI, p_rtriangulation%ndim, NVEMax)
          
          ! Since the discretisations are not uniform, we have to treat each
          ! element individually since it may differ from its predecessor.
          ! However, we may skip the re-initialisation of cubature points, etc.
          ! if the current elements belongs to the same element distribution as 
          ! the last one.

          ! Loop over the patches in the set
          DO ipatch = 1, npatchesInCurrentBlock
            
            ! Loop over elements in patch
            DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
              
              ! Get global element number
              IEL = IelementsInPatch(idx)

              ! Get number of local element distribution
              ilocalElementDistr = p_IelementDistr(IEL)

              ! Check if local element distribution corresponds to the last element 
              ! distribution. Then we don't have to initialise everything again.
              IF (ilocalElementDistr .NE. ilastElementDistr) THEN

                ! Active local element distribution
                p_elementDistribution => p_rdiscrSource%RelementDistribution(ilocalElementDistr)

                ! Get the number of local DOF's for trial functions
                indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)

                ! Get the number of corner vertices of the element
                NVE = elem_igetNVE(p_elementDistribution%itrialElement)

                !---------------------------------------------------------------------
                ! Step 1:  Prepare the source FE space
                !---------------------------------------------------------------------
                
                ! Initialise the cubature formula for the local element distribution
                ! Get cubature weights and point coordinates on the reference element
                CALL cub_getCubPoints(p_elementDistribution%ccubTypeEval, ncubp, Dxi, Domega)

                ! Get from the trial element space the type of coordinate system
                ! that is used in the local element distribution
                icoordSystem = elem_igetCoordSystem(p_elementDistribution%itrialElement)
                               
                ! Check if one of the trial/test elements is nonparametric
                bnonparTrial = elem_isNonparametric(p_elementDistribution%itestElement)
                
                ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                ! p_DcubPtsRef - depending on whether the space is parametric or not.
                IF (bnonparTrial) THEN
                  p_DcubPtsTrial => p_DcubPtsReal
                ELSE
                  p_DcubPtsTrial => p_DcubPtsRef
                END IF
                
                ! Save number of last element distribution
                ilastElementDistr = ilocalElementDistr
              END IF

              ! Put the cubature point coordinates in the right format to the 
              ! cubature-point array.
              ! Initialise all entries in p_DcubPtsRef with the same coordinates -
              ! as the cubature point coordinates are identical on all elements.
              ! Importantly, only those entries are copied which correspond to
              ! the local element type. The remaining entries are left equal to zero.
              DO i=1, ncubp
                DO k=1,SIZE(p_DcubPtsRef,1)
                  ! Could be solved using the TRANSPOSE operator - but often is's 
                  ! faster this way...
                  p_DcubPtsRef(k,i,idx) = Dxi(i,k)
                END DO
              END DO

              ! Update number of sampling points for current patch
              Inpoints(ipatch) = Inpoints(ipatch)+ncubp
              
              ! Store number of corners per element
              IelementNVEInPatch(idx) = NVE

              ! Store number of cubature points per element
              IelementNcubpInPatch(idx) = ncubp

              ! We have the coordinates of the cubature points saved in the
              ! coordinate array from above. Unfortunately for nonparametric
              ! elements, we need the real coordinate.
              ! Furthermore, we anyway need the coordinates of the element
              ! corners and the Jacobian determinants corresponding to
              ! all the points.
              
              ! At first, get the coordinates of the corners of the element.
              CALL trafo_getCoords(&
                  elem_igetTrafoType(p_elementDistribution%itrialElement), &
                  p_rtriangulation, IEL, p_Dcoords(:,:,idx))

              ! Depending on the type of transformation, we must now choose
              ! the mapping between the reference and the real element.
              ! In case we use a nonparametric element as test function, we need the 
              ! coordinates of the points on the real element, too.
              ! Unfortunately, we need the real coordinates of the cubature points
              ! anyway for the function - so calculate them all.
              CALL trafo_calctrafo_mult (p_elementDistribution%ctrafoType, ncubp, &
                  p_Dcoords(:,:,idx), p_DcubPtsRef(:,:,idx), p_Djac(:,:,idx), &
                  p_Ddetj(:,idx), p_DcubPtsReal(:,:,idx))

              !---------------------------------------------------------------------
              ! Step 2:  Perform sampling of consistent gradient values
              !---------------------------------------------------------------------
              
              ! Calculate the derivative of the FE function in the cubature
              ! points: u_h(x,y,z) and save the result to Dcoefficients(:,:,1..3)
              SELECT CASE(p_rtriangulation%ndim)
              CASE (NDIM1D)
                CALL fevl_evaluate_mult (rvectorScalar, p_Dcoords(:,:,idx), &
                    p_Djac(:,:,idx), p_Ddetj(:,idx), &
                    p_elementDistribution%itrialElement, IdofsTrial(:,idx), ncubp, &
                    p_DcubPtsTrial(:,:,idx), DER_DERIV1D_X, Dcoefficients(:,idx,1))
!!$                    DcoefficientsMixed(nnpoints-ncubp+1:nnpoints,1))
                
              CASE (NDIM2D)
                CALL fevl_evaluate_mult (rvectorScalar, p_Dcoords(:,:,idx), &
                    p_Djac(:,:,idx), p_Ddetj(:,idx), &
                    p_elementDistribution%itrialElement, IdofsTrial(:,idx), ncubp, &
                    p_DcubPtsTrial(:,:,idx), DER_DERIV2D_X, Dcoefficients(:,idx,1))
!!$                    DcoefficientsMixed(nnpoints-ncubp+1:nnpoints,1))

                CALL fevl_evaluate_mult (rvectorScalar, p_Dcoords(:,:,idx), &
                    p_Djac(:,:,idx), p_Ddetj(:,idx), &
                    p_elementDistribution%itrialElement, IdofsTrial(:,idx), ncubp, &
                    p_DcubPtsTrial(:,:,idx), DER_DERIV2D_Y, Dcoefficients(:,idx,2))
!!$                    DcoefficientsMixed(nnpoints-ncubp+1:nnpoints,2))
                
              CASE (NDIM3D)
                CALL fevl_evaluate_mult (rvectorScalar, p_Dcoords(:,:,idx), &
                    p_Djac(:,:,idx), p_Ddetj(:,idx), &
                    p_elementDistribution%itrialElement, IdofsTrial(:,idx), ncubp, &
                    p_DcubPtsTrial(:,:,idx), DER_DERIV3D_X, Dcoefficients(:,idx,1))
!!$                    DcoefficientsMixed(nnpoints-ncubp+1:nnpoints,1))
                
                CALL fevl_evaluate_mult (rvectorScalar, p_Dcoords(:,:,idx), &
                    p_Djac(:,:,idx), p_Ddetj(:,idx), &
                    p_elementDistribution%itrialElement, IdofsTrial(:,idx), ncubp, &
                    p_DcubPtsTrial(:,:,idx), DER_DERIV3D_Y, Dcoefficients(:,idx,2))
!!$                    DcoefficientsMixed(nnpoints-ncubp+1:nnpoints,2))
                
                CALL fevl_evaluate_mult (rvectorScalar, p_Dcoords(:,:,idx), &
                    p_Djac(:,:,idx), p_Ddetj(:,idx), &
                    p_elementDistribution%itrialElement, IdofsTrial(:,idx), ncubp, &
                    p_DcubPtsTrial(:,:,idx), DER_DERIV3D_Z, Dcoefficients(:,idx,3))
!!$                    DcoefficientsMixed(nnpoints-ncubp+1:nnpoints,3))
                
              CASE DEFAULT
                CALL output_line('Invalid spatial dimension!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
                CALL sys_halt()
              END SELECT             

            END DO
          END DO   ! End of IPATCH loop
          ilastElementDistr = 0


          !---------------------------------------------------------------------
          ! Step 3: Prepare least-squares fitting
          !---------------------------------------------------------------------

          ! First, get the coordinates of the bounding group for the set of
          ! elements present in each patch. In addition, calculate the coordinates
          ! of the corners of the constant Jacobian patch "elements".
          CALL calc_patchBoundingGroup_mult(IelementsInPatchIdx, npatchesInCurrentBlock, &
              IelementNVEInPatch, p_Dcoords, DpatchBound(:,:,1:npatchesInCurrentBlock))

          ! Reactive current element distribution of the source FE space. This is the 
          ! element distribution of those elements which make up the center of each 
          ! patch. Hence, this element distribution determine the type of element 
          ! used to construct the constant Jacobian "patch" elements.
          p_elementDistribution => p_rdiscrSource%RelementDistribution(icurrentElementDistr)

          ! Initialise the cubature formula for the local element distribution
          ! Get cubature weights and point coordinates on the reference element
          CALL cub_getCubPoints(p_elementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
          
          ! Get the number of local DOF's for trial functions
          indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)

          ! Get the number of corner vertices of the element
          NVE = elem_igetNVE(p_elementDistribution%itrialElement)
          
          ! Get from the trial element space the type of coordinate system
          ! that is used there:
          icoordSystem = elem_igetCoordSystem(p_elementDistribution%itrialElement)

          ! Check if one of the trial/test elements is nonparametric
          bnonparTrial = elem_isNonparametric(p_elementDistribution%itestElement)
          
          ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
          ! p_DcubPtsRef - depending on whether the space is parametric or not.
          IF (bnonparTrial) THEN
            p_DcubPtsTrial => p_DcubPtsReal
          ELSE
            p_DcubPtsTrial => p_DcubPtsRef
          END IF
          
          ! Next, we need to convert the physical coordinates of the curvature points
          ! to the local coordinates of the constant Jacobian "patch" elements
          CALL calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem,&
              npatchesInCurrentBlock, NVE, p_DcubPtsReal, DpatchBound, p_DcubPtsRef)

          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType, &
              nelementsPerBlock, ncubpMax, p_Dcoords, p_DcubPtsRef, p_Djac, p_Ddetj)

          ! Allocate memory for the patch interpolants matrices
          ! Note that the second dimension is DER_FUNC=1 and could be omitted. However,
          ! all routines for simultaneous element evaluation require a 4d-array so that
          ! the array Dpolynomials is artificially created as 4d-array.
          ALLOCATE(Dpolynomials(indofTrial,DER_FUNC,ncubpMax,nelementsPerBlock))

          ! Evaluate the trial functions of the constant Jacobian patch "element" for all
          ! cubature points of the elements present in the patch and store each polynomial 
          ! interpolation in the rectangular patch matrix used for least-squares fitting.
          ! Note that we evaluate over the maximum number of cubature points present in
          ! the patch. Altought some meaningless values may be generated, it is faster to
          ! evaluate all values simultaneously and filter the required data afterwards.
          CALL elem_generic_sim(p_elementDistribution%itrialElement,&
              p_Dcoords, p_Djac, p_Ddetj, BderDest, Dpolynomials, ncubpMax,&
              nelementsPerBlock, p_DcubPtsTrial)

          ! Compute total number of sampling points
          nspoints = SUM(Inpoints(1:npatchesInCurrentBlock))

          ! Allocate memory for the aligned patch data
          ALLOCATE(DcoefficientsMixed(nspoints,p_rtriangulation%ndim))
          ALLOCATE(DpolynomialsMixed(indofTrial*nspoints))
          
          ! Reset total number of sampling/interpolation points
          nspoints = 0; nipoints = 0

          ! Loop over all element indices in the set of patches
          ! and align the polynomials interpolants
          DO idx = 1, IelementsInPatchIdx(npatchesInCurrentBlock+1)-1
            
            ! Get number of cubature points for current element
            ncubp = IelementNcubpInPatch(idx)

            ! Update total number of sampling points
            nspoints = nspoints+ncubp

            ! Apply coefficients for each dimension
            DO idim = 1, p_rtriangulation%ndim
              DcoefficientsMixed(nspoints-ncubp+1:nspoints,idim) =&
                  Dcoefficients(1:ncubp,idx,idim)
            END DO

            ! Loop over all cubature points
            DO icubp = 1, ncubp
              
              ! Update total number of interpolation points
              nipoints = nipoints+indofTrial
              
              ! Align polynomial interpolats
              DpolynomialsMixed(nipoints-indofTrial+1:nipoints) =&
                   Dpolynomials(1:indofTrial,DER_FUNC,icubp,idx)
            END DO
          END DO

          ! Deallocate memory for the unaligned data
          DEALLOCATE(Dpolynomials, Dcoefficients)


          !-----------------------------------------------------------------------
          ! Step 4: Perform least-squares fitting
          !-----------------------------------------------------------------------

          ! Allocate memory for the derivative values
          ALLOCATE(Dderivatives(indofTrial,npatchesInCurrentBlock,p_rtriangulation%ndim))
          
          ! Compute the patch averages by solving $(P^T * P) * x = (P^T) * b$ for x
          CALL calc_patchAverages_mult(IelementsInPatchIdx, npatchesInCurrentBlock, &
              Inpoints, indofTrial, DcoefficientsMixed, DpolynomialsMixed, Dderivatives)

          DEALLOCATE(DpolynomialsMixed, DcoefficientsMixed)


          !-----------------------------------------------------------------------
          ! Step 5: Prepare the destination FE space
          !-----------------------------------------------------------------------

          ! Get local references to domain integration structure
          p_DcubPtsRef =>  rintSubsetDest%p_DcubPtsRef
          p_DcubPtsReal => rintSubsetDest%p_DcubPtsReal
          p_Djac =>        rintSubsetDest%p_Djac
          p_Ddetj =>       rintSubsetDest%p_Ddetj
          p_Dcoords =>     rintSubsetDest%p_DCoords
          
          ! Initialise arrays with zeros
          p_DcubPtsRef  = 0
          p_DcubPtsReal = 0
          p_Djac        = 0
          p_Ddetj       = 0
          p_Dcoords     = 0

          ! Loop over the patches in the set
          DO ipatch = 1, npatchesInCurrentBlock
            
            ! Loop over elements in patch
            DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
              
              ! Get global element number
              IEL = IelementsInPatch(idx)

              ! Get number of local element distribution
              ilocalElementDistr = p_IelementDistr(IEL)

              ! Check if local element distribution corresponds to the last element 
              ! distribution. Then we don't have to initialise everything again.
              IF (ilocalElementDistr .NE. ilastElementDistr) THEN

                ! Active local element distribution
                p_elementDistrDest => p_rdiscrDest%RelementDistribution(ilocalElementDistr)

                ! Get the number of local DOF's for trial functions
                indofDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)

                ! Initialise arrays with zeros
                Dxi    = 0
                Domega = 0

                ! Initialise the cubature formula. That's a special trick here!
                ! In particular, we only need the physical coordinates of the
                ! nodal evaluation points in the source vectors.
                ! For this purpose, we initialise the 'trapezoidal rule' as cubature
                ! formula. Ok, we are not using cubature at all here (thus ignoring
                ! the cubature weights completely), but that way we get the
                ! coordinates of the corners on the reference element automatically --
                ! which coincide with the points where we want to create the gradients!
                !
                ! Note: The returned nlocalDOFsDest will coincide with the number of local DOF's
                ! on each element indofDest!
                CALL calc_cubatureDest(&
                    elem_getPrimaryElement(p_elementDistrDest%itrialElement), &
                    nlocalDOFsDest, Dxi, Domega)
                                
                ! Check if one of the trial/test elements is nonparametric
                bnonparTrial = elem_isNonparametric(p_elementDistrDest%itestElement)
                
                ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                ! p_DcubPtsRef - depending on whether the space is parametric or not.
                IF (bnonparTrial) THEN
                  p_DcubPtsTrial => p_DcubPtsReal
                ELSE
                  p_DcubPtsTrial => p_DcubPtsRef
                END IF
                
                ! Save number of last element distribution
                ilastElementDistr = ilocalElementDistr
              END IF

              ! Put the cubature point coordinates in the right format to the
              ! cubature-point array.
              ! Initialise all entries in p_DcubPtsRef with the same coordinates -
              ! as the cubature point coordinates are identical on all elements
              DO i=1,nlocalDOFsDest
                DO k=1,SIZE(p_DcubPtsRef,1)
                  ! Could be solved using the TRANSPOSE operator - but often is's 
                  ! faster this way...
                  p_DcubPtsRef(k,i,idx) = Dxi(i,k)
                END DO
              END DO
              
              ! We have the coordinates of the cubature points saved in the
              ! coordinate array from above. Unfortunately for nonparametric
              ! elements, we need the real coordinate.
              ! Furthermore, we anyway need the coordinates of the element
              ! corners and the Jacobian determinants corresponding to
              ! all the points.
              
              ! At first, get the coordinates of the corners of the element
              CALL trafo_getCoords (&
                  elem_igetTrafoType(p_elementDistrDest%itrialElement), &
                  p_rtriangulation, IEL , p_Dcoords(:,:,idx))
              
              ! Depending on the type of transformation, we must now choose
              ! the mapping between the reference and the real element.
              ! In case we use a nonparametric element as test function, we need the 
              ! coordinates of the points on the real element, too.
              ! Unfortunately, we need the real coordinates of the cubature points
              ! anyway for the function - so calculate them all.
              CALL trafo_calctrafo_mult (p_elementDistrDest%ctrafoType, nlocalDOFsDest, &
                  p_Dcoords(:,:,idx), p_DcubPtsRef(:,:,idx), p_Djac(:,:,idx), &
                  p_Ddetj(:,idx), p_DcubPtsReal(:,:,idx))
            END DO
          END DO
          ilastElementDistr = 0
          
          ! Next, we need to convert the physical coordinates of the curvature points
          ! to the local coordinates of the constant Jacobian "patch" elements. Note that 
          ! NVE and icoordSystem have not been modified and can be savely used from above.
          CALL calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, npatchesInCurrentBlock, &
              NVE,  p_DcubPtsReal, DpatchBound, p_DcubPtsRef)

          ! We do not need the corner coordinates of the elements in the destination
          ! FE space but that of the "patch" elements in the source FE space
          p_Dcoords => rintSubset%p_Dcoords


          ! It is mandatory that the second dimensions is DER_MAXNDER even though only the
          ! function values (DER_FUNC) are actually needed. During the evaluation of the basis
          ! functions some elements do not check if first derivatives are required. In short,
          ! IF-THEN-ELSE is more costly than just computing some superficial values ;-)
          ! Therefore, we must provide sufficient memory !!!
          ALLOCATE(Dpolynomials(indofTrial,DER_MAXNDER,nlocalDOFsDestMax,nelementsPerBlock))
          Dpolynomials = 0.0_DP

          ! Loop over the patches in the set
          DO ipatch = 1, npatchesInCurrentBlock
            
            ! Loop over elements in patch
            DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
              
              ! Get global element number
              IEL = IelementsInPatch(idx)
              
              ! Get number of local element distribution
              ilocalElementDistr = p_IelementDistr(IEL)
              
              ! Check if local element distribution corresponds to the last element 
              ! distribution. Then we don't have to initialise everything again.
              IF (ilocalElementDistr .NE. ilastElementDistr) THEN
                
                ! Active local element distribution
                p_elementDistrDest => p_rdiscrDest%RelementDistribution(ilocalElementDistr)
                
                ! Get the number of local DOF's for trial functions which coincides with
                nlocalDOFsDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)

                ! Check if one of the trial/test elements is nonparametric
                bnonparTrial = elem_isNonparametric(p_elementDistrDest%itestElement)
                
                ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                ! p_DcubPtsRef - depending on whether the space is parametric or not.
                IF (bnonparTrial) THEN
                  p_DcubPtsTrial => p_DcubPtsReal
                ELSE
                  p_DcubPtsTrial => p_DcubPtsRef
                END IF

                ! Save number of last element distribution
                ilastElementDistr = ilocalElementDistr
              END IF
                
              ! Calculate the transformation from the reference element to the real one
              CALL trafo_calctrafo_mult(p_elementDistrDest%ctrafoType, nlocalDOFsDest, &
                  p_Dcoords(:,:,idx), p_DcubPtsRef(:,:,idx), p_Djac(:,:,idx), p_Ddetj(:,idx))

              ! Evaluate the basis functions for the cubature points of the destination FE space
              CALL elem_generic_mult(p_elementDistribution%itrialElement, &
                  p_Dcoords(:,:,idx), p_Djac(:,:,idx), p_Ddetj(:,idx), &
                  BderDest, Dpolynomials(:,:,:,idx), nlocalDOFsDest, &
                  p_DcubPtsTrial(:,:,idx))
            END DO
          END DO
          ilastElementDistr = 0

          !---------------------------------------------------------------------
          ! Step 6: Evaluate the averaged derivative values at the cubature
          !         points of the destination FE space and scatter them to
          !         the global degrees of freedom
          !---------------------------------------------------------------------

          SELECT CASE (p_rtriangulation%ndim)
          CASE(NDIM1D)
            
            ! Loop over the patches in the set
            DO ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1

                ! Get global element number
                IEL = IelementsInPatch(idx)
                
                ! Get number of local element distribution
                ilocalElementDistr = p_IelementDistr(IEL)
                
                ! Check if local element distribution corresponds to the last element 
                ! distribution. Then we don't have to initialise everything again.
                IF (ilocalElementDistr .NE. ilastElementDistr) THEN
                  
                  ! Active local element distribution
                  p_elementDistrDest => p_rdiscrDest%RelementDistribution(ilocalElementDistr)
                  
                  ! Get the number of local DOF's for trial functions which coincides with
                  nlocalDOFsDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)

                  ! Check if one of the trial/test elements is nonparametric
                  bnonparTrial = elem_isNonparametric(p_elementDistrDest%itestElement)
                  
                  ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                  ! p_DcubPtsRef - depending on whether the space is parametric or not.
                  IF (bnonparTrial) THEN
                    p_DcubPtsTrial => p_DcubPtsReal
                  ELSE
                    p_DcubPtsTrial => p_DcubPtsRef
                  END IF

                  ! Save number of last element distribution
                  ilastElementDistr = ilocalElementDistr
                END IF

                ! Loop over local degrees of freedom
                DO ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  DO j = 1, indofTrial
                    Dval(1) = Dval(1) + Dderivatives(j,ipatch,1) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  END DO
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                END DO
              END DO
            END DO
            ilastElementDistr = 0

          CASE(NDIM2D)
            
            ! Loop over the patches in the set
            DO ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1

                ! Get global element number
                IEL = IelementsInPatch(idx)
                
                ! Get number of local element distribution
                ilocalElementDistr = p_IelementDistr(IEL)
                
                ! Check if local element distribution corresponds to the last element 
                ! distribution. Then we don't have to initialise everything again.
                IF (ilocalElementDistr .NE. ilastElementDistr) THEN
                  
                  ! Active local element distribution
                  p_elementDistrDest => p_rdiscrDest%RelementDistribution(ilocalElementDistr)
                  
                  ! Get the number of local DOF's for trial functions which coincides with
                  nlocalDOFsDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)
                  
                  ! Save number of last element distribution
                  ilastElementDistr = ilocalElementDistr
                END IF

                ! Loop over local degrees of freedom
                DO ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  DO j = 1, indofTrial
                    Dval(1:2) = Dval(1:2) + Dderivatives(j,ipatch,1:2) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  END DO
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                END DO
              END DO
            END DO
            ilastElementDistr = 0

          CASE(NDIM3D)
            
            ! Loop over the patches in the set
            DO ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              DO idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1

                ! Get global element number
                IEL = IelementsInPatch(idx)
                
                ! Get number of local element distribution
                ilocalElementDistr = p_IelementDistr(IEL)
                
                ! Check if local element distribution corresponds to the last element 
                ! distribution. Then we don't have to initialise everything again.
                IF (ilocalElementDistr .NE. ilastElementDistr) THEN
                  
                  ! Active local element distribution
                  p_elementDistrDest => p_rdiscrDest%RelementDistribution(ilocalElementDistr)
                  
                  ! Get the number of local DOF's for trial functions which coincides with
                  nlocalDOFsDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)
                  
                  ! Save number of last element distribution
                  ilastElementDistr = ilocalElementDistr
                END IF

                ! Loop over local degrees of freedom
                DO ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  DO j = 1, indofTrial
                    Dval(1:3) = Dval(1:3) + Dderivatives(j,ipatch,1:3) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  END DO
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  p_DzDeriv(IdofsDest(ipoint,idx)) = p_DzDeriv(IdofsDest(ipoint,idx)) + Dval(3)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                END DO
              END DO
            END DO
            ilastElementDistr = 0

          CASE DEFAULT
            CALL output_line('Invalid spatial dimension!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
            CALL sys_halt()
          END SELECT

          ! Release memory
          CALL domint_doneIntegration(rintSubset)
          CALL domint_doneIntegration(rintSubsetDest)

          ! Deallocate temporary memory
          DEALLOCATE(Dderivatives)
          DEALLOCATE(Dpolynomials)
          DEALLOCATE(IdofsDest)
          DEALLOCATE(IdofsTrial)
          DEALLOCATE(IelementNVEInPatch)
          DEALLOCATE(IelementNcubpInPatch)
        END IF
        
        ! Deallocate temporary memory
        DEALLOCATE(IelementsInPatch)

      END DO   ! End of PATCHset loop
      
      ! Deallocate temporary memory
      IF (.NOT.bisUniform) DEALLOCATE(Inpoints)
      DEALLOCATE(DpatchBound)
      DEALLOCATE(IelementsInPatchIdx)

    END DO   ! End of icurrentElementDistr loop

    ! We are nearly done. The final thing: divide the calculated derivatives by the
    ! number of elements adjacent to each vertex. That closes the calculation
    ! of the 'mean' of the derivatives.
    DO i=1,SIZE(p_DxDeriv)
      ! Div/0 should not occur, otherwise the triangulation is crap as there's a point
      ! not connected to any element!
      p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
      p_DyDeriv(i) = p_DyDeriv(i) / p_IcontributionsAtDOF(i)
    END DO

    ! Free temporary memory
    CALL storage_free(h_IcontributionsAtDOF)
    
  CONTAINS
    
    ! Here, some auxiliary working routines follow

    !**************************************************************
    ! Calculate the "bounding group" for a set of patches
    !
    ! Each patch consists of multiple elements which are adjacent to 
    ! each other so that they cover some connected part of the domain.
    ! This routine determines the physical coordinates of the constant
    ! Jacobian patch "element" that has its local axes parallel to the
    ! global axes and completely surrounds the elements of the patch.
    ! Moreover, this routine converts the physical coordinates of the
    ! evaluation points in the destination space to the local coordinates
    ! in the constant Jacobian patch "element".

    SUBROUTINE calc_patchBoundingGroup_mult(IelementsInPatchIdx, &
        npatches, IelemNVE, DpointsReal, DpointsBound)
      
      ! Index vector and list of elements in patch
      INTEGER, DIMENSION(:), INTENT(IN)         :: IelementsInPatchIdx

      ! Number of elements in patch
      INTEGER, INTENT(IN)                       :: npatches

      ! Array with numbers of corners per element
      INTEGER, DIMENSION(:), INTENT(IN)         :: IelemNVE

      ! Physical coordinates of the corner nodes of all elements in the patch
      REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: DpointsReal

      ! Physical coordinates of the bounds of all "patch" elements
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT)   :: DpointsBound

      ! local variables
      INTEGER :: ipatch,idx,idxFirst,idxLast
      REAL(DP) :: xmin,ymin,xmax,ymax

      SELECT CASE (SIZE(DpointsReal,1))

      CASE (NDIM1D)
        
        ! Loop over all patches
        DO ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1
          
          ! Determine minimum/maximum coordinates of the patch
          xmin = MINVAL(DpointsReal(1,:,idxFirst:idxLast))
          xmax = MAXVAL(DpointsReal(1,:,idxFirst:idxLast))
          
          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(1,2,ipatch) = xmax
          
          ! Store physical coordinates of patch element corners
          DpointsReal(1,1,idxFirst:idxLast) = xmin
          DpointsReal(1,2,idxFirst:idxLast) = xmax
        END DO
        
      CASE (NDIM2D)

        ! Loop over all patches
        DO ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1
          
          ! Determine minimum/maximum coordinates of the patch
          xmin = MINVAL(DpointsReal(1,:,idxFirst:idxLast))
          xmax = MAXVAL(DpointsReal(1,:,idxFirst:idxLast))
          ymin = MINVAL(DpointsReal(2,:,idxFirst:idxLast))
          ymax = MAXVAL(DpointsReal(2,:,idxFirst:idxLast))
          
          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(2,1,ipatch) = ymin
          DpointsBound(1,2,ipatch) = xmax
          DpointsBound(2,2,ipatch) = ymax
          
          ! Loop over all elements
          DO idx = idxFirst, idxLast
            SELECT CASE(IelemNVE(idx))
            CASE(TRIA_NVETRI2D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(1,2,idx) = xmax+ymax-ymin
              DpointsReal(2,2,idx) = ymin
              DpointsReal(1,3,idx) = xmin
              DpointsReal(2,3,idx) = xmax+ymax-xmin

            CASE (TRIA_NVEQUAD2D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(1,2,idx) = xmax
              DpointsReal(2,2,idx) = ymin
              DpointsReal(1,3,idx) = xmax
              DpointsReal(2,3,idx) = ymax
              DpointsReal(1,4,idx) = xmin
              DpointsReal(2,4,idx) = ymax

            CASE DEFAULT
              CALL output_line ('Invalid number of vertices per elements!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_mult')
              CALL sys_halt()
            END SELECT
          END DO
        END DO
        
      CASE DEFAULT
        CALL output_line ('Invalid number of spatial dimensions!', &
            OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_mult')
        CALL sys_halt()
      END SELECT
    END SUBROUTINE calc_patchBoundingGroup_mult
    
    !**************************************************************
    ! Calculate the "bounding group" for a set of patches
    !
    ! Each patch consists of multiple elements which are adjacent to 
    ! each other so that they cover some connected part of the domain.
    ! This routine determines the physical coordinates of the constant
    ! Jacobian patch "element" that has its local axes parallel to the
    ! global axes and completely surrounds the elements of the patch.
    ! Moreover, this routine converts the physical coordinates of the
    ! evaluation points in the destination space to the local coordinates
    ! in the constant Jacobian patch "element".

    SUBROUTINE calc_patchBoundingGroup_sim(IelementsInPatchIdx, &
        npatches, NVE, DpointsReal, DpointsBound)
      
      ! Index vector and list of elements in patch
      INTEGER, DIMENSION(:), INTENT(IN)         :: IelementsInPatchIdx

      ! Number of elements in patch
      INTEGER, INTENT(IN)                       :: npatches

      ! Number of vertices per element
      INTEGER, INTENT(IN)                       :: NVE

      ! Physical coordinates of the corner nodes of all elements in the patch
      REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: DpointsReal

      ! Physical coordinates of the bounds of all "patch" elements
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT)   :: DpointsBound

      ! local variables
      INTEGER :: ipatch,idxFirst,idxLast
      REAL(DP) :: xmin,ymin,xmax,ymax

      SELECT CASE (NVE)

      CASE (TRIA_NVELINE1D)

        ! Simple: Just find minimal/maximal value
        
        ! Loop over all patches
        DO ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Determine minimum/maximum coordinates of the patch
          xmin = MINVAL(DpointsReal(1,:,idxFirst:idxLast))
          xmax = MAXVAL(DpointsReal(1,:,idxFirst:idxLast))

          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(1,2,ipatch) = xmax

          ! Store physical coordinates of patch element corners
          DpointsReal(1,1,idxFirst:idxLast) = xmin
          DpointsReal(1,2,idxFirst:idxLast) = xmax
        END DO

      CASE(TRIA_NVETRI2D)
        
        ! Tricky: Find minimal/maximal value and compute the lower-right
        ! and upper-left corner by hand.

        ! Loop over all patches
        DO ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Determine minimum/maximum coordinates of the patch
          xmin = MINVAL(DpointsReal(1,:,idxFirst:idxLast))
          xmax = MAXVAL(DpointsReal(1,:,idxFirst:idxLast))
          ymin = MINVAL(DpointsReal(2,:,idxFirst:idxLast))
          ymax = MAXVAL(DpointsReal(2,:,idxFirst:idxLast))
          
          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(2,1,ipatch) = ymin
          DpointsBound(1,2,ipatch) = xmax
          DpointsBound(2,2,ipatch) = ymax

          ! Store physical coordinates of patch element corners
          DpointsReal(1,1,idxFirst:idxLast) = xmin
          DpointsReal(2,1,idxFirst:idxLast) = ymin
          DpointsReal(1,2,idxFirst:idxLast) = xmax+ymax-ymin
          DpointsReal(2,2,idxFirst:idxLast) = ymin
          DpointsReal(1,3,idxFirst:idxLast) = xmin
          DpointsReal(2,3,idxFirst:idxLast) = xmax+ymax-xmin
        END DO
        
      CASE (TRIA_NVEQUAD2D)
        
        ! Simple: Just find minimal/maximal value

        ! Loop over all patches
        DO ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Determine minimum/maximum coordinates of the patch
          xmin = MINVAL(DpointsReal(1,:,idxFirst:idxLast))
          xmax = MAXVAL(DpointsReal(1,:,idxFirst:idxLast))
          ymin = MINVAL(DpointsReal(2,:,idxFirst:idxLast))
          ymax = MAXVAL(DpointsReal(2,:,idxFirst:idxLast))
                  
          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(2,1,ipatch) = ymin
          DpointsBound(1,2,ipatch) = xmax
          DpointsBound(2,2,ipatch) = ymax

          ! Store physical coordinates of patch element corners
          DpointsReal(1,1,idxFirst:idxLast) = xmin
          DpointsReal(2,1,idxFirst:idxLast) = ymin
          DpointsReal(1,2,idxFirst:idxLast) = xmax
          DpointsReal(2,2,idxFirst:idxLast) = ymin
          DpointsReal(1,3,idxFirst:idxLast) = xmax
          DpointsReal(2,3,idxFirst:idxLast) = ymax
          DpointsReal(1,4,idxFirst:idxLast) = xmin
          DpointsReal(2,4,idxFirst:idxLast) = ymax
        END DO
        
      CASE DEFAULT
        CALL output_line ('Invalid number of vertices per elements!', &
            OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_sim')
        CALL sys_halt()
      END SELECT
    END SUBROUTINE calc_patchBoundingGroup_sim

    !**************************************************************
    ! Transform physical coordinates to local coordinates of the 
    ! constant Jacobian "patch" elements which is uniquely 
    ! determined by its minimum/maximum values.
    ! 
    ! Note that for quadrilateral/hexahedral elements the x-, y- and
    ! z-coordinates are stored whereas for triangular/tetrahedral
    ! elements barycentric coordinates are adopted.
    
    SUBROUTINE calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, &
        npatches, NVE, DpointsReal, DpointsBound, DpointsRef)

      ! Index vector and list of elements in patch
      INTEGER, DIMENSION(:), INTENT(IN)       :: IelementsInPatchIdx

      ! Coordinate system identifier. One of the TRAFO_CS_xxxx constants. Defines
      ! the type of the coordinate system that is used for specifying the coordinates
      ! on the reference element.
      INTEGER, INTENT(IN)                     :: icoordSystem
      
      ! Number of elements in patch
      INTEGER, INTENT(IN)                     :: npatches

      ! Number of vertices per element
      INTEGER, INTENT(IN)                     :: NVE
      
      ! Coordinates of the points on the physical element.
      REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: DpointsReal

      ! Coordinates of the bounding group of each patch
      REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: DpointsBound

      ! Coordinates of the points on the reference elements.
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: DpointsRef
      
      ! local variables
      INTEGER :: ipatch,idxFirst,idxLast
      REAL(DP) :: xmin,ymin,xmax,ymax,daux

      ! Which coordinate system should be applied
      SELECT CASE (icoordSystem)

      CASE (TRAFO_CS_BARY2DTRI)
        ! This coordinate system can only be applied if NVE=TRIA_NVETRI2D
        IF (NVE .NE. TRIA_NVETRI2D) THEN
          CALL output_line('Invalid number of corner vertices per element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'calc_localTrafo_sim')
          CALL sys_halt()
        END IF
        
        ! Loop over all patches
        DO ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Get minimum/maximum coordinates of the patch
          xmin = DpointsBound(1,1,ipatch)
          ymin = DpointsBound(2,1,ipatch)
          xmax = DpointsBound(1,2,ipatch)
          ymax = DpointsBound(2,2,ipatch)

          daux = xmax-xmin+ymax-ymin

          ! Transform physical coordinates to local ones by means
          ! of a unit coordinate transformation
          DpointsRef(3,:,idxFirst:idxLast) = (DpointsReal(2,:,idxFirst:idxLast)-ymin)/daux
          DpointsRef(2,:,idxFirst:idxLast) = (DpointsReal(1,:,idxFirst:idxLast)-xmin)/daux
          DpointsRef(1,:,idxFirst:idxLast) = 1._DP-DpointsRef(2,:,idxFirst:idxLast) &
                                                  -DpointsRef(3,:,idxFirst:idxLast)
        END DO
        
      CASE (TRAFO_CS_REF2DTRI,TRAFO_CS_REF2DQUAD,&
            TRAFO_CS_REAL2DTRI,TRAFO_CS_REAL2DQUAD,TRAFO_CS_REF1D)
        
        ! How many corner vertices do we have?
        SELECT CASE (NVE)
        CASE (TRIA_NVELINE1D)
          
          ! Loop over all patches
          DO ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Get minimum/maximum coordinates of the patch
            xmin = DpointsBound(1,1,ipatch)
            xmax = DpointsBound(1,2,ipatch)
            
            ! Transform physical coordinates to local ones by means
            ! of a natural coordinate transformation
            DpointsRef(1,:,idxFirst:idxLast) = &
                (2.0_DP*DpointsReal(1,:,idxFirst:idxLast)-(xmax+xmin))/(xmax-xmin)
          END DO

        CASE(TRIA_NVETRI2D)
          
          ! Loop over all patches
          DO ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Get minimum/maximum coordinates of the patch
            xmin = DpointsBound(1,1,ipatch)
            ymin = DpointsBound(2,1,ipatch)
            xmax = DpointsBound(1,2,ipatch)
            ymax = DpointsBound(2,2,ipatch)
            
            ! Transform physical coordinates to local ones by means
            ! of a unit coordinate transformation
            DpointsRef(1,:,idxFirst:idxLast) = &
                (DpointsReal(1,:,idxFirst:idxLast)-(xmax+ymax-ymin))/(xmax-xmin+ymax-ymin)
            DpointsRef(2,:,idxFirst:idxLast) = &
                (DpointsReal(2,:,idxFirst:idxLast)-(xmax+ymax-xmin))/(xmax-xmin+ymax-ymin)
          END DO
          
        CASE (TRIA_NVEQUAD2D)
          
          ! Loop over all patches
          DO ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Get minimum/maximum coordinates of the patch
            xmin = DpointsBound(1,1,ipatch)
            ymin = DpointsBound(2,1,ipatch)
            xmax = DpointsBound(1,2,ipatch)
            ymax = DpointsBound(2,2,ipatch)
            
            ! Transform physical coordinates to local ones by means
            ! of a natural coordinate transformation
            DpointsRef(1,:,idxFirst:idxLast) = &
                (2.0_DP*DpointsReal(1,:,idxFirst:idxLast)-(xmax+xmin))/(xmax-xmin)
            DpointsRef(2,:,idxFirst:idxLast) = &
                (2.0_DP*DpointsReal(2,:,idxFirst:idxLast)-(ymax+ymin))/(ymax-ymin)
          END DO
          
        CASE DEFAULT
          CALL output_line ('Invalid number of vertices per elements!', &
              OU_CLASS_ERROR,OU_MODE_STD,'calc_localTrafo_sim')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unknown coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'calc_localTrafo_sim')
        CALL sys_halt()
      END SELECT
    END SUBROUTINE calc_localTrafo_sim
    
    !**************************************************************    
    ! Calculate the averaged gradient values
    ! In principal, it suffices to solve the linear system
    ! $$ (P^T * P) * x = (P^T) * b$$
    ! for the unknown $x$. However, least-squares fitting is known
    ! to be ill-conditions. An alternative is suggested in
    !
    ! J.E. Akin, Finite Element Analysis with Error Estimation
    !
    ! Akin suggests to perform a singular value decomposition (SVD)
    ! of the matrix (P^T * P) and perform back substitution.
    ! That's exactly, what is done in this subroutine.
    
    SUBROUTINE calc_patchAverages_mult(IelementsInPatchIdx, npatches, &
        Inpoints, indofTrial, Dcoefficients, Dpolynomials, Dderivatives)
      
      ! Index vector and list of elements in patch
      INTEGER, DIMENSION(:), INTENT(IN)                   :: IelementsInPatchIdx

      ! Number of elements in patch
      INTEGER, INTENT(IN)                                 :: npatches

      ! Array with numbers of sampling points per patch
      INTEGER, DIMENSION(:), INTENT(IN)                   :: Inpoints

      ! Number of local degrees of freedom for test functions
      INTEGER, INTENT(IN)                                 :: indofTrial

      ! Vector with consistent gradient values
      REAL(DP), DIMENSION(:,:), INTENT(IN)                :: Dcoefficients

      ! Vector with aligned polynomial interpolants
      REAL(DP), DIMENSION(:), INTENT(INOUT)               :: Dpolynomials

      ! Smoothe gradient values
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT)             :: Dderivatives

      ! local variables
      INTEGER  :: i,ipatch,icoeffFirst,icoeffLast,ipolyFirst,ipolyLast
      REAL(DP), DIMENSION(indofTrial,indofTrial)     :: Dv
      REAL(DP), DIMENSION(indofTrial)                :: Dd
      

      ! Initialise position index
      icoeffFirst = 1
      
      ! Loop over all patches
      DO ipatch = 1, npatches
        
        ! Compute absolute position of coefficient and polynomials
        icoeffLast = icoeffFirst+Inpoints(ipatch)-1
        ipolyFirst = indofTrial*(icoeffFirst-1)+1
        ipolyLast  = indofTrial*icoeffLast

        ! Compute factorisation for the singular value decomposition
        CALL mprim_SVD_factorise(Dpolynomials(ipolyFirst:ipolyLast),&
            indofTrial, Inpoints(ipatch), Dd, Dv, .TRUE.)
        
        ! Perform back substitution for all componenets
        DO i = 1, SIZE(Dderivatives,3)
          CALL mprim_SVD_backsubst1(Dpolynomials(ipolyFirst:ipolyLast), &
              indofTrial, Inpoints(ipatch), Dd, Dv, Dderivatives(:,ipatch,i), &
              Dcoefficients(icoeffFirst:icoeffLast,i), .TRUE.)
        END DO

        ! Update position index
        icoeffFirst = icoeffFirst+Inpoints(ipatch)
      END DO
    END SUBROUTINE calc_patchAverages_mult

    !**************************************************************    
    ! Calculate the averaged gradient values
    ! In principal, it suffices to solve the linear system
    ! $$ (P^T * P) * x = (P^T) * b$$
    ! for the unknown $x$. However, least-squares fitting is known
    ! to be ill-conditions. An alternative is suggested in
    !
    ! J.E. Akin, Finite Element Analysis with Error Estimation
    !
    ! Akin suggests to perform a singular value decomposition (SVD)
    ! of the matrix (P^T * P) and perform back substitution.
    ! That's exactly, what is done in this subroutine.
    
    SUBROUTINE calc_patchAverages_sim(IelementsInPatchIdx, npatches, &
        ncubp, indofTrial, Dcoefficients, Dpolynomials, Dderivatives)

      ! Index vector and list of elements in patch
      INTEGER, DIMENSION(:), INTENT(IN)              :: IelementsInPatchIdx

      ! Number of elements in patch
      INTEGER, INTENT(IN)                            :: npatches

      ! Number of cubature points
      INTEGER, INTENT(IN)                            :: ncubp

      ! Number of local degrees of freedom for test functions
      INTEGER, INTENT(IN)                            :: indofTrial

      ! Vector with consistent gradient values
      !
      !   Dcoefficients(ncubp, nelemPerBlock, ndim)
      !
      REAL(DP), DIMENSION(:,:,:), INTENT(IN)         :: Dcoefficients

      ! Rectangular matrix with polynomial interpolants
      !
      !   Dpolynomials(indofTrial, 1, ncubp, nelemPerBlock)
      !
      REAL(DP), DIMENSION(:,:,:,:), INTENT(INOUT)    :: Dpolynomials

      ! Smoothed gradient values
      !
      !   Derivatives(indofTrial, npatches, ndim)
      !
      REAL(DP), DIMENSION(:,:,:), INTENT(OUT)        :: Dderivatives

      ! local variables
      INTEGER  :: i,ipatch,idxFirst,idxLast,npoints

      REAL(DP), DIMENSION(indofTrial,indofTrial)     :: Dv
      REAL(DP), DIMENSION(indofTrial)                :: Dd
      
      ! Loop over all patches
      DO ipatch = 1, npatches
        
        ! Get first and last element index of patch
        idxFirst = IelementsInPatchIdx(ipatch)
        idxLast  = IelementsInPatchIdx(ipatch+1)-1
        npoints  = ncubp*(idxLast-idxFirst+1)

        ! Compute factorisation for the singular value decomposition
        CALL mprim_SVD_factorise(Dpolynomials(:,:,:,idxFirst:idxLast),&
            indofTrial, npoints, Dd, Dv, .TRUE.)
        
        ! Perform back substitution for all componenets
        DO i = 1, SIZE(Dderivatives,3)
          CALL mprim_SVD_backsubst2(Dpolynomials(:,:,:,idxFirst:idxLast), &
              indofTrial, npoints, Dd, Dv, Dderivatives(:,ipatch,i), &
              Dcoefficients(:,idxFirst:idxLast,i), .TRUE.)
        END DO
      END DO
    END SUBROUTINE calc_patchAverages_sim

    !**************************************************************    
    ! Initialise the cubature formula for the destination FE space

    SUBROUTINE calc_cubatureDest(ieltyp,ncubp, Dxi, Domega)
      
      ! Element type identifier.
      INTEGER, INTENT(IN) :: ieltyp

      ! Number of cubature points
      INTEGER, INTENT(OUT) :: ncubp

      ! Coordinates of the cubature points
      REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dxi
      
      ! Cubature weights of the cubature points
      REAL(DP), DIMENSION(:), INTENT(OUT) :: Domega
      
      SELECT CASE (ieltyp)
      CASE (EL_P0)
        CALL cub_getCubPoints(CUB_G1_T, ncubp, Dxi, Domega)
      CASE (EL_Q0)
        CALL cub_getCubPoints(CUB_G1X1, ncubp, Dxi, Domega)
      CASE (EL_P1)
        CALL cub_getCubPoints(CUB_TRZ_T, ncubp, Dxi, Domega)
      CASE (EL_Q1)
        CALL cub_getCubPoints(CUB_TRZ, ncubp, Dxi, Domega)
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
        
        ncubp = 6
        
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
        
        ncubp = 9
        
      CASE DEFAULT 
        CALL output_line ('Unsupported FE space in destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'calc_cubatureDest')
        CALL sys_halt()
      END SELECT
    END SUBROUTINE calc_cubatureDest
  END SUBROUTINE ppgrd_calcGradSuperPatchRecov

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppgrd_calcGradLimAvgP1Q1cnf (rvectorScalar,rvectorGradient)

!<description>
    ! Calculates the recovered gradient of a scalar finite element function
    ! by means of the limited gradient averaging technique suggested by 
    ! M. Möller and D. Kuzmin. Supports conformal discretisations in arbitrary
    ! spatial dimensions with $P_1$ and $Q_1$ finite elements mixed in the
    ! source and destination vectors.
!</description>

!<input>
    ! The FE solution vector. Represents a scalar FE function.
    TYPE(t_vectorScalar), INTENT(IN)         :: rvectorScalar
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
    LOGICAL :: bisQ1, bisP1, bisQ1T, bisP1T, bisDifferent
    LOGICAL :: bnonparTrial
    INTEGER(I32) :: IELmax, IELset, idof

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
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: IdofsTrial
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: IdofsDest
    
    ! Pointers to the X-, Y- and Z-derivative vector
    REAL(DP), DIMENSION(:), POINTER :: p_DxDeriv, p_DyDeriv, p_DzDeriv
    
    ! Number of elements in the current element distribution
    INTEGER(PREC_ELEMENTIDX) :: NEL

    ! Pointer to an array that counts if an edge has been visited.
    INTEGER :: h_IcontributionsAtDOF
    INTEGER(I32), DIMENSION(:), POINTER :: p_IcontributionsAtDOF
    
    ! Discretisation structures for the source- and destination vector(s)
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrSource, p_rdiscrDest

    ! Dimension of triangulation must be less or equal than number of subvectors
    ! in the block vector.
    
    IF (.NOT. ASSOCIATED(rvectorScalar%p_rspatialDiscretisation)) THEN
      CALL output_line ('No discretisation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF
    
    IF ((rvectorScalar%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND.&
        (rvectorScalar%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
      CALL output_line ('Only uniform and conformal discretisations supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF
    
    IF (.NOT. ASSOCIATED(rvectorScalar%p_rspatialDiscretisation%p_rtriangulation)) THEN
      CALL output_line ('No triangulation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF
    
    IF (rvectorScalar%p_rspatialDiscretisation%p_rtriangulation%ndim .GT. &
        rvectorGradient%nblocks) THEN
      CALL output_line ('Dimension of destination vector not large enough!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF
    
    ! There must be given discretisation structures in the destination vector.
    IF (.NOT. ASSOCIATED(rvectorScalar%p_rspatialDiscretisation)) THEN
      CALL output_line ('No discretisation attached to the destination vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF

    ! The source vector must be either pure Q1 or pure P1 or mixed Q1/P1
    bisQ1 = .FALSE.
    bisP1 = .FALSE.
    bisDifferent = .FALSE.

    p_rdiscrSource => rvectorScalar%p_rspatialDiscretisation
    DO j=1,p_rdiscrSource%inumFESpaces
      SELECT CASE (elem_getPrimaryElement (p_rdiscrSource%RelementDistribution(j)%itrialElement))
      CASE (EL_Q1)
        bisQ1 = .TRUE.
      CASE (EL_P1)
        bisP1 = .TRUE.
      CASE DEFAULT
        bisDifferent = .TRUE.
      END SELECT
    END DO

    IF (bisDifferent) THEN
      CALL output_line ('Only Q1, and P1 supported as&
          & discretisation for the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF

    ! The destination vector must be either pure Q1T or pure P1T or mixed Q1T/P1T
    bisQ1T = .FALSE.
    bisP1T = .FALSE.
    bisDifferent = .FALSE.
    
    DO i=1,MIN(rvectorGradient%nblocks,&
        rvectorScalar%p_rspatialDiscretisation%p_rtriangulation%ndim,rvectorGradient%nblocks)
      p_rdiscrDest => rvectorGradient%p_rblockDiscretisation%RspatialDiscretisation(i)
      DO j=1,p_rdiscrDest%inumFESpaces
        SELECT CASE (elem_getPrimaryElement (p_rdiscrDest%RelementDistribution(j)%itrialElement))
        CASE (EL_Q1T)
          bisQ1T = .TRUE.
        CASE (EL_P1T)
          bisP1T = .TRUE.
        CASE DEFAULT
          bisDifferent = .TRUE.
        END SELECT
      END DO
    END DO

    IF (bisDifferent) THEN
      CALL output_line ('Only Q1T, and P1T supported as&
          & discretisation for the destination vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END IF

    ! Get the discretisation structures of the source- and destination space.
    ! Note that we assume here that the X- and Y-derivative is discretised
    ! the same way!
    p_rdiscrSource => rvectorScalar%p_rspatialDiscretisation
    p_rdiscrDest   => rvectorGradient%p_rblockDiscretisation%RspatialDiscretisation(1)
    
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
                     
    ! Get pointers to the derivative destination vector.
    SELECT CASE(p_rtriangulation%ndim)
    CASE (NDIM1D)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      CALL lalg_clearVectorDble (p_DxDeriv)

    CASE (NDIM2D)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      CALL lalg_clearVectorDble (p_DxDeriv)
      CALL lalg_clearVectorDble (p_DyDeriv)

    CASE (NDIM3D)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      CALL lsyssc_getbase_double (rvectorGradient%RvectorBlock(3),p_DzDeriv)
      CALL lalg_clearVectorDble (p_DxDeriv)
      CALL lalg_clearVectorDble (p_DyDeriv)
      CALL lalg_clearVectorDble (p_DzDeriv)

    CASE DEFAULT
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END SELECT

    ! Array that allows the calculation about the number of elements
    ! meeting in a vertex, based onm DOF's.
    CALL storage_new ('ppgrd_calcGradLimAvgP1Q1cnf','p_IcontributionsAtDOF',&
                      dof_igetNDofGlob(p_rdiscrDest),&
                      ST_INT, h_IcontributionsAtDOF, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int (h_IcontributionsAtDOF,p_IcontributionsAtDOF)
    
    ! Evaluate the first derivative of the FE functions.

    Bder = .FALSE.
    SELECT CASE (p_rtriangulation%ndim)
    CASE (NDIM1D)
      Bder(DER_DERIV1D_X) = .TRUE.

    CASE (NDIM2D)
      Bder(DER_DERIV2D_X) = .TRUE.
      Bder(DER_DERIV2D_Y) = .TRUE.

    CASE (NDIM3D)
      Bder(DER_DERIV3D_X) = .TRUE.
      Bder(DER_DERIV3D_Y) = .TRUE.
      Bder(DER_DERIV3D_Z) = .TRUE.

    CASE DEFAULT
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      CALL sys_halt()
    END SELECT

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    DO icurrentElementDistr = 1,p_rdiscrSource%inumFESpaces
    
      ! Activate the current element distribution
      p_elementDistribution => p_rdiscrSource%RelementDistribution(icurrentElementDistr)
      p_elementDistrDest => p_rdiscrDest%RelementDistribution(icurrentElementDistr)
    
      ! If the element distribution is empty, skip it
      IF (p_elementDistribution%NEL .EQ. 0) CYCLE
    
      ! Get the number of local DOF's for trial functions
      ! in the source and destination vector.
      indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
      indofDest = elem_igetNDofLoc(p_elementDistrDest%itrialElement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_elementDistribution%itrialElement)

      ! Initialise the cubature formula.
      ! That's a special trick here! The FE space of the destination vector
      ! is either P1 or Q1. We create the gradients in the midpoints of the edges
      ! by taking the limited average of the gradients of the source vector!

      SELECT CASE (elem_getPrimaryElement(p_elementDistrDest%itrialElement))
      CASE (EL_P1T)
        CALL cub_getCubPoints(CUB_G3_T, nlocalDOFsDest, Dxi, Domega)

      CASE (EL_Q1T)
        Dxi(1,1) =  0.0_DP
        Dxi(1,2) = -1.0_DP
        Dxi(2,1) =  1.0_DP
        Dxi(2,2) =  0.0_DP
        Dxi(3,1) =  0.0_DP
        Dxi(3,2) =  1.0_DP
        Dxi(4,1) = -1.0_DP
        Dxi(4,2) =  0.0_DP
        
        nlocalDOFsDest = 4

      CASE DEFAULT
        CALL output_line ('Unsupported FE space in destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
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

      ! Destination space is either P1 or Q1.
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

      ! Get the number of elements in the element distribution.
      NEL = p_elementDistribution%NEL
      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                                p_IelementList)

      ! Loop over the elements - blockwise.
      DO IELset = 1, NEL, PPGRD_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = MIN(NEL,IELset-1+PPGRD_NELEMSIM)
      
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
        
        CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_elementDistribution%itrialElement, IdofsTrial, &
              nlocalDOFsDest, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV_X,&
              Dderivatives(:,1:IELmax-IELset+1_I32,1))        

        CALL fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_elementDistribution%itrialElement, IdofsTrial, &
              nlocalDOFsDest, INT(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV_Y,&
              Dderivatives(:,1:IELmax-IELset+1_I32,2))

        ! Sum up the derivative values in the destination vector.
        ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
        ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
        ! local DOF's in the destination space -- in that order!

        SELECT CASE (p_rtriangulation%ndim)
        CASE (NDIM1D)
          
          DO i=1,IELmax-IELset+1
            DO j=1,nlocalDOFsDest

              idof = IdofsDest(j,i)

              IF (p_IcontributionsAtDOF(idof) .EQ. 0) THEN 
                p_DxDeriv(idof) = Dderivatives(j,i,1)
                p_IcontributionsAtDOF(idof) = 1
              ELSE
                p_DxDeriv(idof) = (SIGN(0.5_DP,p_DxDeriv(idof)) + &
                                   SIGN(0.5_DP,Dderivatives(j,i,1))) * &
                                   MIN(0.5_DP*ABS(p_DxDeriv(idof)+Dderivatives(j,i,1)), &
                                       2._DP*ABS(p_DxDeriv(idof)), &
                                       2._DP*ABS(Dderivatives(j,i,1)))
              END IF
            END DO
          END DO

        CASE (NDIM2D)

          DO i=1,IELmax-IELset+1
            DO j=1,nlocalDOFsDest

              idof = IdofsDest(j,i)

              IF (p_IcontributionsAtDOF(idof) .EQ. 0) THEN 
                p_DxDeriv(idof) = Dderivatives(j,i,1)
                p_DyDeriv(idof) = Dderivatives(j,i,2)
                p_IcontributionsAtDOF(idof) = 1
              ELSE
                p_DxDeriv(idof) = (SIGN(0.5_DP,p_DxDeriv(idof)) + &
                                   SIGN(0.5_DP,Dderivatives(j,i,1))) * &
                                   MIN(0.5_DP*ABS(p_DxDeriv(idof)+Dderivatives(j,i,1)), &
                                       2._DP*ABS(p_DxDeriv(idof)), &
                                       2._DP*ABS(Dderivatives(j,i,1)))
                p_DyDeriv(idof) = (SIGN(0.5_DP,p_DyDeriv(idof)) + &
                                   SIGN(0.5_DP,Dderivatives(j,i,2))) * &
                                   MIN(0.5_DP*ABS(p_DyDeriv(idof)+Dderivatives(j,i,2)), &
                                       2._DP*ABS(p_DyDeriv(idof)), &
                                       2._DP*ABS(Dderivatives(j,i,2)))
              END IF
            END DO
          END DO

        CASE (NDIM3D)

          DO i=1,IELmax-IELset+1
            DO j=1,nlocalDOFsDest

              idof = IdofsDest(j,i)

              IF (p_IcontributionsAtDOF(idof) .EQ. 0) THEN 
                p_DxDeriv(idof) = Dderivatives(j,i,1)
                p_DyDeriv(idof) = Dderivatives(j,i,2)
                p_DzDeriv(idof) = Dderivatives(j,i,3)
                p_IcontributionsAtDOF(idof) = 1
              ELSE
                p_DxDeriv(idof) = (SIGN(0.5_DP,p_DxDeriv(idof)) + &
                                   SIGN(0.5_DP,Dderivatives(j,i,1))) * &
                                   MIN(0.5_DP*ABS(p_DxDeriv(idof)+Dderivatives(j,i,1)), &
                                       2._DP*ABS(p_DxDeriv(idof)), &
                                       2._DP*ABS(Dderivatives(j,i,1)))
                p_DyDeriv(idof) = (SIGN(0.5_DP,p_DyDeriv(idof)) + &
                                   SIGN(0.5_DP,Dderivatives(j,i,2))) * &
                                   MIN(0.5_DP*ABS(p_DyDeriv(idof)+Dderivatives(j,i,2)), &
                                       2._DP*ABS(p_DyDeriv(idof)), &
                                       2._DP*ABS(Dderivatives(j,i,2)))
                p_DzDeriv(idof) = (SIGN(0.5_DP,p_DzDeriv(idof)) + &
                                   SIGN(0.5_DP,Dderivatives(j,i,3))) * &
                                   MIN(0.5_DP*ABS(p_DzDeriv(idof)+Dderivatives(j,i,3)), &
                                       2._DP*ABS(p_DzDeriv(idof)), &
                                       2._DP*ABS(Dderivatives(j,i,3)))
              END IF
            END DO
          END DO

        CASE DEFAULT
          CALL output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
          CALL sys_halt()
        END SELECT
        
      END DO ! IELset
      
      ! Release memory
      CALL domint_doneIntegration(rintSubset)
      
      DEALLOCATE(Dderivatives)
      DEALLOCATE(IdofsDest)
      DEALLOCATE(IdofsTrial)
      
    END DO ! icurrentElementDistr

    ! Release temp data
    CALL storage_free (h_IcontributionsAtDOF)

  END SUBROUTINE ppgrd_calcGradLimAvgP1Q1cnf

END MODULE
