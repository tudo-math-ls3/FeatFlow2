!##############################################################################
!# ****************************************************************************
!# <name> jumpstabilisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines that implement the jump stabilisation.
!#
!# The following routines can be found in this module:
!#
!# 1.) jstab_calcUEOJumpStabilisation
!#     -> Modifies a scalar matrix to add the unified edge-oriented
!#        jump stabilisation term.
!# 
!# Some auxiliary routines:
!#
!# 1.) jstab_ueoJumpStabil2d_m_unidble
!#     -> The actual computation routine for 2D domains.
!#
!# </purpose>
!##############################################################################

MODULE jumpstabilisation

  USE fsystem
  USE genoutput
  USE linearsystemscalar
  USE linearsystemblock
  USE cubature
  USE domainintegration
  USE bilinearformevaluation
  USE element
  USE elementpreprocessing
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE jstab_calcUEOJumpStabilisation (&
      rmatrix,dgamma,dgammastar,dtheta,ccubType,dnu)

!<description>
  ! Edge oriented stabilisation technique. 
!</description>

!<input>
  ! Stabilisation parameter. Standard=0.01
  REAL(DP), INTENT(IN) :: dgamma
  
  ! 2nd stabilisation parametr. Standard=dgamma=0.01
  REAL(DP), INTENT(IN) :: dgammastar
  
  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  REAL(DP), INTENT(IN) :: dtheta
  
  ! 1D cubature formula to use for line integration.
  ! Standard = CUB_G2_1D.
  INTEGER, INTENT(IN) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  REAL(DP), INTENT(IN) :: dnu
!</input>

!<inputoutput>
  ! Scalar matrix to be modified. The stabilisation term is added to rmatrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format 
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    IF ((rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
        (rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
      CALL output_line ('Unsupported matrix format.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      CALL sys_halt()
    END IF

    IF (rmatrix%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
      CALL output_line ('Unsupported discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      CALL sys_halt()
    END IF

    ! Modify the matrix.
    CALL jstab_ueoJumpStabil2d_m_unidble (&
        rmatrix,dgamma,dgammastar,dtheta,ccubType,dnu)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE jstab_ueoJumpStabil2d_m_unidble ( &
      rmatrixScalar,dgamma,dgammastar,dtheta,ccubType,dnu)
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! $$< Ju,v > = \sum_E \max(\gamma^{*}*\nu*h_E, \gamma h_E^2) 
  !            \int_E [grad u] [grad v] ds$$
  !
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
  !
  ! WARNING: For edge oriented stabilisation, the underlying matrix rmatrix
  !   must have an extended stencil! The matrix structure must be set up with
  !   the BILF_MATC_EDGEBASED switch!!!
!</description>
  
!<input>
  ! Stabilisation parameter. Standard=0.01
  REAL(DP), INTENT(IN) :: dgamma
  
  ! 2nd stabilisation parametr. Standard=dgamma=0.01
  REAL(DP), INTENT(IN) :: dgammastar
  
  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  REAL(DP), INTENT(IN) :: dtheta
  
  ! 1D cubature formula to use for line integration
  ! Standard = CUB_G2_1D.
  INTEGER, INTENT(IN) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  REAL(DP), INTENT(IN) :: dnu
!</input>
  
!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(PREC_DOFIDX) :: irow, jcol, idof
  INTEGER(PREC_EDGEIDX) :: IMT
  INTEGER(PREC_VERTEXIDX) :: ivt1,ivt2,NVT
  INTEGER(PREC_ELEMENTIDX) :: IEL
  LOGICAL :: bIdenticalTrialAndTest
  INTEGER :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  REAL(DP) :: dval,dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTestTempl
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialTempl
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2), TARGET :: IdofsTest
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2), TARGET :: IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IdofsTrial
  
  ! Arrays saving the local DOF numbers belonging to the global
  ! DOF numbers in IdofsTest/IdofsTrial  
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2),TARGET :: IlocalDofsTest,IlocalDofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IlocalDofsTrial
  
  ! Renumbering strategy for local DOF's
  INTEGER, DIMENSION(EL_MAXNBAS), TARGET :: IlocalDofTestRenum,IlocalDofTrialRenum
  INTEGER, DIMENSION(:), POINTER :: p_IlocalDofTrialRenum
  
  ! Number of local DOF's on the patch
  INTEGER :: ndofTest,ndofTrial
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrialPerElement, indofTestPerElement
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! Some triangulation arrays we need frequently
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IelementsAtEdge
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_relementDistribution
  
  ! Underlying discretisation structure
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation

  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj

  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Dentry

  ! Type of transformation from the reference to the real element 
  INTEGER :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  INTEGER(I32) :: cevaluationTag
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  TYPE(t_domainIntSubset) :: rintSubset
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

  ! Cubature point weights
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp,icubp
  
  ! Derivative specifiers
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrial, BderTest
  
  ! Whther the viscosity ís constant or not
  LOGICAL :: bconstViscosity
  
    ! Currently we support only constant viscosity
    bconstViscosity = .TRUE.

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar%p_rspatialDiscretisation%p_rtriangulation
    p_rdiscretisation => rmatrixScalar%p_rspatialDiscretisation
    
    ! Get Kvert, Kadj, Kmid, Kmel,...
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
                                p_IneighboursAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IelementsAtEdge,&
                                p_IelementsAtEdge)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)
    CALL storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)
                                
    NVT = p_rtriangulation%NVT
    
    ! Get KLD, KCol...
    CALL lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    CALL lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    CALL lsyssc_getbase_double (rmatrixScalar,p_Da)
    
    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistribution(1)

    ! Get the number of local DOF's for trial and test functions
    indofTrialPerElement = elem_igetNDofLoc(p_relementDistribution%itrialElement)
    indofTestPerElement = elem_igetNDofLoc(p_relementDistribution%itestElement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%itrialElement)
    
    ! Get the number of corner vertices of the element
    IF (NVE .NE. elem_igetNVE(p_relementDistribution%itestElement)) THEN
      CALL output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_unidble')
      CALL sys_halt()
    END IF
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference line [-1,1]
    CALL cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega)

    ! Allocate arrays accepting cubature point coordinates.
    ! We have ncubp cubature points on every egde.
    ! DcubPtsRef saves all possible combinations, which edge of
    ! one element might interact with one edge of one neighbour
    ! element.
    ALLOCATE(DcubPtsRefOnAllEdges(NDIM2D,ncubp,NVE))
    
    ! Put the 1D cubature points from to all of the edges on the
    ! reference element, so we can access them quickly.
    DO iedge = 1,NVE
      CALL trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
      DO i=1,ncubp
        DcubPtsRefOnAllEdges(1,i,iedge) = Dxi2D(i,1)
        DcubPtsRefOnAllEdges(2,i,iedge) = Dxi2D(i,2)
      END DO
    END DO
    
    ! We have always 2 elements in each patch...
    IELcount = 2
      
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%itrialElement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP,IELcount))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    ALLOCATE(Kentry(indofTestPerElement*2*indofTrialPerElement*2))
    ALLOCATE(Dentry(indofTestPerElement*2*indofTrialPerElement*2))
    
    ! Allocate memory for obtaining DOF's:
    ALLOCATE(IdofsTrialTempl(indofTrialPerElement,2))
    ALLOCATE(IdofsTestTempl(indofTestPerElement,2))
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    !
    ! We allocate space for 3 instead of 2 elements. The reason is that
    ! we later permute the values of the basis functions to get
    ! a local numbering on the patch. That's also the reason, we allocate
    ! not indofTestPerElement elements, but even indofTestPerElement,
    ! which is more than enough space to hold the values of the DOF's of
    ! a whole element patch.
    
    ALLOCATE(DbasTest(indofTestPerElement*2, &
            elem_getMaxDerivative(p_relementDistribution%itestElement),&
            ncubp,3))
    ALLOCATE(DbasTrial(indofTrialPerElement*2,&
            elem_getMaxDerivative(p_relementDistribution%itrialElement), &
            ncubp,3))
             
    ! Test if trial/test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = &
      p_relementDistribution%itrialElement .EQ. p_relementDistribution%itestElement
      
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    ! The same for p_IlocalDofsTrial and the permutation array p_IlocalDofTrialRenum.
    ! Let p_DbasTrial point either to DbasTrial or DbasTest, depending on
    ! whether the spaces are identical.
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
      p_IlocalDofsTrial => IlocalDofsTest
      p_IlocalDofTrialRenum => IlocalDofTestRenum
      p_DbasTrial => DbasTest
    ELSE
      p_IdofsTrial => IdofsTrial
      p_IlocalDofsTrial => IlocalDofsTrial
      p_IlocalDofTrialRenum => IlocalDofTrialRenum
      p_DbasTrial => DbasTrial
    END IF
    
    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%itrialElement)
    cevaluationTag = IOR(cevaluationTag,&
                    elem_getEvaluationTag(p_relementDistribution%itestElement))

    ! Don't calculate coordinates on the reference element -- we do this manually.                    
    cevaluationTag = IAND(cevaluationTag,EL_EVLTAG_REFPOINTS)


    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    BderTrial = .FALSE.
    BderTrial(DER_DERIV_X) = .TRUE.
    BderTrial(DER_DERIV_Y) = .TRUE.

    BderTest = .FALSE.
    BderTest(DER_DERIV_X) = .TRUE.
    BderTest(DER_DERIV_Y) = .TRUE.

    ! We loop through all edges
    DO IMT = 1,p_rtriangulation%NMT
    
      ! Check if we have 1 or 2 elements on the current edge
      IF (p_IelementsAtEdge (2,IMT) .EQ. 0) THEN
      
        ! This only happens if we have a boundary edge.
        ! The boundary handling is... doing nothing! So we can skip
        ! the handling of that edge completely!
        CYCLE
      
      END IF
      
      ! Fill the basis function arrays with 0. Essential, as only parts
      ! of them are overwritten later.
      DbasTest = 0.0_DP
      DbasTrial = 0.0_DP
      
      ! On an example, we now show the relationship between all the DOF's
      ! we have to consider in the current situation. We have two elements,
      ! let's say IEL1 and IEL2, with their local and global DOF's
      ! (example for Q1~):
      !
      !   local DOF's on the           corresponding global
      !   element and on the patch     DOF's
      !
      !    +----4----+----7----+       +----20---+----50---+
      !    |    4    |    4    |       |         |         |
      !    |         |         |       |         |         |
      !    1 1     3 3 1     3 6       10        60        70
      !    |         |         |       |         |         |
      !    |    2    |    2    |       |         |         |
      !    +----2----+----3----+       +----40---+----30---+
      !        IEL1      IEL2              IEL1      IEL2
      !
      !
      ! On every element, we have 4 local DOF's (1..4). On the other hand,
      ! we have "local DOF's on the patch" (1..7), ehich we call "patch DOF's"
      ! from now on. To every local DOF, there belongs a global DOF
      ! (10,20,30,... or whatever), which gives the coefficient of the basis
      ! function.
      !
      ! Numbering that way, the local DOF's of IEL1 obviously coincide 
      ! with the first couple of local DOF's of the element patch! Only
      ! the local DOF's of element IEL2 make trouble, as they have another
      ! numbering.
      
      ! Get the global DOF's of the 1 or two elements
      CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  .TRUE., IdofsTestTempl)
                                   
      ! Some of the DOF's on element 2 may coincide with DOF's on element 1.
      ! More precisely, some must coincide! Therefore, we now have to collect the
      ! DOF's uniquely and to figure out, which local DOF's of element 2
      ! must renumbered to the appropriate local patch-DOF (like in the
      ! example, where local DOF 1 of element IEL2 must be renumbered
      ! to local patch DOF 3!
      !
      ! As the first couple of local DOF's of IEL1 coincide with the local
      ! DOF's of the patch, we can simply copy them:
      
      ndofTest = indofTestPerElement
      IdofsTest(1:ndofTest) = IdofsTestTempl(1:ndofTest,1)
      
      ! Furthermore, we read IdofsTestTempl and store the DOF's there in IdofsTest,
      ! skipping all DOF's we already have and setting up the renumbering strategy
      ! of local DOF's on element IEL2 to patch DOF's.
      
      skiploop: DO IDOFE = 1,indofTestPerElement
        
        ! Do we have the DOF? 
        idof = IdofsTestTempl(IDOFE,IELcount)
        DO JDOFE=1,ndofTest
          IF (IdofsTest(JDOFE) .EQ. idof) THEN
            ! Yes, we have it.
            ! That means, the local DOF idof has to be mapped to the
            ! local patch dof...
            IlocalDofTestRenum (IDOFE) = JDOFE
            
            ! Proceed with the next one.
            CYCLE skiploop
          END IF
        END DO
        
        ! We don't have that DOF! Append it to IdofsTest.
        ndofTest = ndofTest+1
        IdofsTest(ndofTest) = idof
        IlocalDofTestRenum (IDOFE) = ndofTest
        
        ! Save also the number of the local DOF.
        ! Note that the global DOF's in IdofsTestTempl(1..indofTestPerElement)
        ! belong to the local DOF's 1..indofTestPerElement -- in that order!
        IlocalDofsTest(ndofTest) = IDOFE
        
      END DO skiploop
        
      ! If the trial functions are different, get those DOF's
      IF (.NOT. bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementsAtEdge (1:IELcount,IMT), &
                                    .FALSE., IdofsTrialTempl)

        ! Collect the DOF's uniquely in IdofsTrial.
        ! The DOF's of the first element are taken as they are.
        
        ndofTrial = indofTrialPerElement
        IdofsTrial(1:ndofTrial) = IdofsTrialTempl(1:ndofTrial,1)
        
        ! From the 2nd element on, we have to check which
        ! DOF's already appear in the first element.
        
        skiploop2: DO IDOFE = 1,indofTrialPerElement
          
          ! Do we have the DOF? 
          idof = IdofsTrialTempl(IDOFE,IELcount)
          DO JDOFE=1,ndofTrial
            IF (IdofsTrial(JDOFE) .EQ. idof) THEN
              ! Yes, we have it.
              ! That means, the local DOF idof has to be mapped to the
              ! local patch dof...
              IlocalDofTrialRenum (IDOFE) = JDOFE
              
              ! Proceed with the next one.
              CYCLE skiploop2
            END IF
          END DO
          
          ! We don't have that DOF! Append it to IdofsTrial.
          ndofTrial = ndofTrial+1
          IdofsTrial(ndofTrial) = idof
          IlocalDofTrialRenum (IDOFE) = ndofTrial
          
          ! Save also the number of the local DOF.
          ! Note that the global DOF's in IdofsTrialTempl(1..indofTrialPerElement)
          ! belong to the local DOF's 1..indofTrialPerElement -- in that order!
          IlocalDofsTrial(ndofTrial) = IDOFE
          
        END DO skiploop2
          
      ELSE
       
        ! Number of DOF's in test and trial space the same, as DOF's are the same.
        ndofTrial = ndofTest
        
      END IF
      
      ! Now we know: Our 'local' matrix (consisting of only these DOF's we just
      ! calculated) is a ndofsTest*ndofsTrial matrix.
      !
      ! Now extract the corresponding entries from the matrix.
      ! Kentry is an index where the entries of the local matrix Dentry
      ! (which we build up later) can be found in the global matrix.
      !
      ! The test space gives us the rows; loop through them
      
      DO IDOFE = 0,ndofTest-1
      
        ! The 'local' DOF IDOFE corresponds in the global matrix to line...
        irow = IdofsTest (1+IDOFE)
        
        ! Loop through that line to find the columns, indexed by the local
        ! DOF's in the trial space.
        trialspaceloop: DO JDOFE = 1,ndofTrial
          
          DO jcol = p_Kld(irow),p_Kld(irow+1)-1
            
            IF (p_Kcol(jcol) .EQ. p_IdofsTrial(JDOFE)) THEN
            
              ! Found! Put as destination pointer to Kentry
              Kentry (IDOFE*ndofTrial+JDOFE) = jcol
              
              ! Clear local matrix
              Dentry (IDOFE*ndofTrial+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              CYCLE trialspaceloop
            
            END IF
            
          END DO

          CALL output_line ('Matrix invalid! Trial-DOF not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_unidble')
          CALL sys_halt()
            
        END DO trialspaceloop
      
      END DO ! JDOFE
      
      ! Now we can set up the local matrix in Dentry. Later, we'll plug it into
      ! the global matrix using the positions in Kentry.
      !
      ! The next step is to evaluate the basis functions in the cubature
      ! points on the edge. To compute the jump, this has to be done
      ! for both elements on the edge. 
      !
      ! Figure out which edge on the current element is IMT.
      ! We need this local numbering later to determine on which edge we
      ! have to place the cubature points.
      DO i = 1,IELcount
        
        IEL = p_IelementsAtEdge(i,IMT)
        DO iedge = 1,NVE
          IF (p_IedgesAtElement (iedge,IEL) .EQ. IMT+NVT) EXIT
        END DO
        
        ! Copy the coordinates of the corresponding cubature points
        ! to DcubPtsEval. We calculated the cubature points on the
        ! reference element in advance, so we can simply take them.
        
        p_DcubPtsRef(:,:,i) = DcubPtsRefOnAllEdges (:,:,iedge)
        
      END DO

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      CALL elprep_prepareSetForEvaluation (rintSubset%revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementsAtEdge (1:IELcount,IMT), &
          ctrafoType,DpointsRef=p_DcubPtsRef)
      p_Ddetj => rintSubset%revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      CALL elem_generic_sim2 (p_relementDistribution%itestElement, &
          rintSubset%revalElementSet, BderTest, DbasTest)

      ! Apply the permutation of the local DOF's on the test functions
      ! on element 2. The numbers of the local DOF's on element 1
      ! coincides with the numbers of the local DOF's on the patch.
      ! Those on element 2, we have to renumber according to the permutation
      ! so that they are in the correct order according to the DOF's on the patch.
      !
      ! We copy the values of the basis functions to the space in
      ! p_DcubPtsTest which is reserved for the 3rd element!
      ! That way, DbasTest has:
      !
      ! local DOF on element   :       1       2       3       4
      !
      ! Global DOF on element 1:      10      40      60      20
      ! DbasTest(1..4,*,*,1):    d1psi10 d1psi40 d1psi60 d1psi20
      !
      ! Global DOF on elemenr 2:      60      30      70      50
      ! DbasTest(1..4,*,*,2):    d2psi60 d2psi30 d2psi70 d2psi50
      ! 
      ! and will be transformed to:
      !
      ! local patch-DOF:            1       2       3       4       5      6        7
      ! global DOF:                10      40      60      20      30     70       50
      ! DbasTest(1..7,*,*,1): d1psi10 d1psi40 d1psi60 d1psi20     0.0    0.0      0.0
      ! DbasTest(1..7,*,*,2): --------------------------unused-----------------------
      ! DbasTest(1..7,*,*,3):     0.0     0.0 d2psi60     0.0 d2psi30 d2psi70 d2psi50
      !
      ! ("d1psi10" = grad(psi_10) on element 1, 
      !  "psi_10" = test function on global DOF #10)
      !
      ! That space DbasTest(:,:,:,3) is unused up to now, initialised by 0.0 and 
      ! now overwritten.
      DbasTest (IlocalDofTestRenum(1:indofTestPerElement),:,:,3) &
        = DbasTest (1:indofTestPerElement,:,:,2)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      IF (.NOT. bidenticalTrialAndTest) THEN
        CALL elem_generic_sim2 (p_relementDistribution%itrialElement, &
            rintSubset%revalElementSet, BderTrial, DbasTrial)

        ! Apply the renumbering for the 2nd element, store the result in the
        ! space of the 3rd element.          
        DbasTrial (IlocalDofTrialRenum(1:indofTrialPerElement),:,:,3) &
          = DbasTrial (1:indofTrialPerElement,:,:,2)
            
      END IF
      
      ! What do we have now? When looking at the example, we have:
      !
      ! ndofTest = 7
      ! IdofsTest(1..7)             = global DOF's in local numbering 1..7
      ! DbasTest(1..7,*,1..ncubp,1) = values of basis functions on element 1
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF's only appearing at element 2
      ! DbasTest(1..7,*,1..ncubp,3) = values of basis functions on element 2
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF's only appearing at element 1
      !
      ! Now we can start to integrate using this.
      
      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit inverval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        SQRT ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 &
            + (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2 )
      dedgeweight = 0.5_DP * dedgelength
      
      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E max(gammastar*nu*h_E, gamma*h_E^2) int_E [grad u] [grad v] ds
      dcoeff = MAX(dgammastar * dnu * dedgelength, &
                   dgamma * dedgelength**2)
      
      ! Now we have the values of the basis functions in all the cubature 
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      DO IDOFE = 1,ndofTest
      
        ! Loop through the trial basis functions
        DO JDOFE = 1,ndofTrial
    
          dval = 0.0_DP
          
          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          DO icubp = 1,ncubp
          
            ! [ grad phi ]   ( jump in the derivative of trial basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dphidx = p_DbasTrial (JDOFE,DER_DERIV_X,icubp,1) &
                  - p_DbasTrial (JDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dphidy = p_DbasTrial (JDOFE,DER_DERIV_Y,icubp,1) &
                  - p_DbasTrial (JDOFE,DER_DERIV_Y,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsidx = DbasTest (IDOFE,DER_DERIV_X,icubp,1) &
                  - DbasTest (IDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dpsidy = DbasTest (IDOFE,DER_DERIV_Y,icubp,1) &
                  - DbasTest (IDOFE,DER_DERIV_Y,ncubp-icubp+1,3)
            
          
            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dval = dval + Domega(icubp) * dedgeweight * &
                          (dphidx*dpsidx + dphidy*dpsidy)
          
          END DO

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndofTrial+JDOFE) = &
            Dentry ((IDOFE-1)*ndofTrial+JDOFE) + dcoeff*dval

        END DO
      
      END DO
      
      ! Incorporate our "local" system matrix
      ! into the global matrix. The position of each entry DENTRY(X,Y)    
      ! in the global matrix array A was saved in element Kentry(X,Y)
      ! before.                                                      
      ! Kentry gives the position of the additive contributions in Dentry.
      ! The entry is weighted by the current dtheta, which is usually
      ! the weighting parameter of the corresponding THETA-scheme of a
      ! nonstationary simulation. For stationary simulations, dtheta is typically
      ! 1.0 which includes the local matrix into the global one directly.)
      
      DO IDOFE = 0,ndofTest-1
        DO JDOFE = 1,ndofTrial
        
          p_Da(Kentry(IDOFE*ndofTrial+JDOFE)) = &
            p_Da(Kentry(IDOFE*ndofTrial+JDOFE)) + &
            dtheta * Dentry (IDOFE*ndofTrial+JDOFE)
        
        END DO
      END DO
      
      ! Proceed with next edge
    
    END DO ! IMT

    ! Clean up allocated arrays and memory.
    DEALLOCATE(DcubPtsRefOnAllEdges)

    CALL elprep_releaseElementSet(rintSubset%revalElementSet)
    DEALLOCATE(p_DcubPtsRef)
    
    DEALLOCATE(DbasTest)
    DEALLOCATE(DbasTrial)

    DEALLOCATE(Kentry) 
    DEALLOCATE(Dentry)

    DEALLOCATE(IdofsTestTempl) 
    DEALLOCATE(IdofsTrialTempl)

  END SUBROUTINE

END MODULE
